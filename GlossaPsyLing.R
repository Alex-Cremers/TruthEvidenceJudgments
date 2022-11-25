library(tidyverse)
library(rstan) # For Bayesian data analysis
library(bridgesampling) # To compute marginal log-likelihood and Bayes factors
library(bayestestR) # For HDI credible intervals
library(xtable) # For latex tables
library(multidplyr) # for parallel dplyr (not crucial, only used to speed up generation of predictions for graphs)

# Make use of all cores:
options(mc.cores = parallel::detectCores())

# Show full tibbles
options(tibble.width=Inf)

# Set up cluster for mutlidplyr:
cluster <- new_cluster(parallel::detectCores()-1)
cluster_library(cluster, "tidyverse")
cluster_copy(cluster,c("ci"))


# Standard error function:
se <- function(x) sd(x)/sqrt(length(x))

# Function to make prompt names prettier:
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

# Load the data
Data <- read_csv("AnonData.csv")

# Error rate on each control question:
Data %>%
  summarize(
    Comprehension1=1-mean(Comprehension1),
    Comprehension2=1-mean(Comprehension2),
    Comprehension3=1-mean(Comprehension3)
  )

# Mean error rate:
1-summarise(Data,mean(Comprehension1+Comprehension2+Comprehension3,na.rm=T)/3)[[1]]

# How many participants passed how many controls?
table(with(Data,Comprehension1+Comprehension2+Comprehension3))

# Keep only participants who passed all three controls
Data <- Data %>%
  filter(Comprehension1+Comprehension2+Comprehension3==3) %>%
  mutate(
    NormAnswer = case_when(
      AnswerType=="Binary" ~ Answer,
      AnswerType=="Likert" ~ (Answer-1)/6
    ),
    Evidence = case_when(
      AnnaTotal == 2 ~ 0,
      AnnaTotal == 4 ~ 0.125,
      AnnaTotal == 7 ~ 1/2,
      AnnaTotal ==10 ~ 0.875,
      AnnaTotal ==12 ~ 1
    ),
    Truth = (Winner=="ANNA")
  )



# Clean up prompt factor:
Data$Prompt <- factor(Data$Prompt,levels=c("TRUE","CORRECT","RIGHT","TRUSTWORTHY","NATURAL","JUSTIFIED"))
Data$Prompt <- factor(Data$Prompt,labels=capwords(levels(Data$Prompt),strict=T))

# How many observations per condition? (we'll display this in the graph below)
SampleSizes <- Data %>% 
  group_by(AnswerType,Prompt) %>%
  summarize(N=n())

Data %>% 
  filter(AnswerType=="Binary") %>%
  group_by(Prompt,Truth) %>%
  summarize(N=n()) %>%
  pull(N) %>%
  chisq.test()
Data %>% 
  filter(AnswerType=="Likert") %>%
  group_by(Prompt,Truth) %>%
  summarize(N=n()) %>%
  pull(N) %>%
  chisq.test()


###############
# First graphs
###############

# Get the levels in the right order:
levels_binary <- paste0(SampleSizes$Prompt[SampleSizes$AnswerType=="Binary"]," (N=",SampleSizes$N[SampleSizes$AnswerType=="Binary"],")")
levels_Likert <- paste0(SampleSizes$Prompt[SampleSizes$AnswerType=="Likert"]," (N=",SampleSizes$N[SampleSizes$AnswerType=="Likert"],")")

# Binary responses (Figure 2a):
pdf("Graphs/Figure2a.pdf",width=10, height=6)
Data %>%
  group_by(Prompt,AnswerType,Truth,Evidence) %>%
  summarize(
    mean=mean(NormAnswer),
    se=se(NormAnswer)
  ) %>%
  ungroup() %>%
  filter(AnswerType=="Binary") %>%
  left_join(SampleSizes) %>%
  mutate(Prompt = factor(paste0(Prompt," (N=",N,")"),levels=levels_binary)) %>%
  ggplot(aes(x=Evidence,y=mean,ymin=mean-se,ymax=mean+se,color=Truth,fill=Truth))+
  facet_wrap(.~Prompt)+
  geom_point()+
  geom_errorbar()+
  geom_line()+
  theme_bw()+
  ylab("Mean answer (SE)")+
  coord_cartesian(ylim=c(-0.01,1.01))+
  scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
  scale_color_manual(values=c("#F5793A","#0F2080"),name="Truth:",labels=c("False","True"),aesthetics = c("colour", "fill"),guide=guide_legend(keywidth=2,reverse=T))
dev.off()

# Likert responses (Figure 2b):
pdf("Graphs/Figure2b.pdf",width=10, height=6)
Data %>%
  group_by(Prompt,AnswerType,Truth,Evidence) %>%
  summarize(
    mean=mean(Answer),
    se=se(Answer)
  ) %>%
  ungroup() %>%
  filter(AnswerType=="Likert") %>%
  left_join(SampleSizes) %>%
  mutate(Prompt = factor(paste0(Prompt," (N=",N,")"),levels=levels_Likert)) %>%
  ggplot(aes(x=Evidence,y=mean,ymin=mean-se,ymax=mean+se,color=Truth,fill=Truth))+
  facet_wrap(.~Prompt)+
  geom_point()+
  geom_errorbar()+
  geom_line()+
  theme_bw()+
  coord_cartesian(ylim=c(1-0.01,7.01))+
  scale_y_continuous(name="Mean answer (SE)",breaks = 1:7)+
  scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
  scale_color_manual(values=c("#F5793A","#0F2080"),name="Truth:",labels=c("False","True"),aesthetics = c("colour", "fill"),guide=guide_legend(keywidth=2,reverse=T))
dev.off()

############################
# Some more data processing
############################

# Define centered predictors for modelling:
Data <- Data %>%
  mutate(
    scaled_evidence = as.numeric(scale(Evidence)),
    centered_truth = ifelse(Winner=="ANNA",0.5,-0.5)
  )

StatDataBinary <- Data %>%
  filter(AnswerType=="Binary") %>%
  mutate(Answer = as.logical(Answer))

StatDataLikert <- Data %>%
  filter(AnswerType=="Likert") %>%
  mutate(Answer = factor(Answer))


#############################
# Illustrating non-linearity
#############################

# Reviewers and editor point out that there seems to be non-linear effects each ends of the Evidence scale
# Here is an illustration of how adding a fixed offset at each end captures the non-linearities:

# Without offset:
model <- Data %>%
  filter(AnswerType=="Binary"&Prompt=="Justified") %>%
  glm(Answer~Evidence*Truth,data=.,family="binomial")
linear_prediction <- expand_grid(Truth=c(F,T),Evidence=c(0,seq(.125,.875,by=.005),1)) %>%
  filter(!(Evidence==0&Truth)&!(Evidence==1&!Truth)) %>%
  mutate(y=plogis(predict(model,newdata=.)))


Data %>%
  filter(AnswerType=="Binary"&Prompt=="Justified") %>%
  group_by(Evidence,Truth) %>%
  summarize(mean_answer=mean(Answer),se_answer=se(Answer)) %>%
  ggplot(aes(x=Evidence,y=mean_answer,ymin=mean_answer-se_answer,ymax=mean_answer+se_answer,color=Truth))+
  geom_point()+
  geom_errorbar()+
  geom_line(data=linear_prediction,aes(x=Evidence,y=y,color=Truth),size=1.2,inherit.aes = F)+
  theme_bw()+
  coord_cartesian(ylim=c(-0.01,1.01))+
  scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"),name="Truth:",labels=c("False","True"),aesthetics = c("colour", "fill"),guide=guide_legend(keywidth=2,reverse=T))

# With offset:
model <- Data %>%
  filter(AnswerType=="Binary"&Prompt=="Justified") %>%
  glm(Answer~Evidence*Truth+I(Evidence==0)+I(Evidence==1),data=.,family="binomial")
non_linear_prediction <- expand_grid(Truth=c(F,T),Evidence=c(0,seq(.05,.95,by=.005),1)) %>%
  filter(!(Evidence==0&Truth)&!(Evidence==1&!Truth)) %>%
  mutate(y=plogis(predict(model,newdata=.)))

Data %>%
  filter(AnswerType=="Binary"&Prompt=="Justified") %>%
  group_by(Evidence,Truth) %>%
  summarize(mean_answer=mean(Answer),se_answer=se(Answer)) %>%
  ggplot(aes(x=Evidence,y=mean_answer,ymin=mean_answer-se_answer,ymax=mean_answer+se_answer,color=Truth))+
  geom_point()+
  geom_errorbar()+
  geom_line(data=non_linear_prediction,aes(x=Evidence,y=y,color=Truth),size=1.2,inherit.aes = F)+
  theme_bw()+
  coord_cartesian(ylim=c(-0.01,1.01))+
  scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"),name="Truth:",labels=c("False","True"),aesthetics = c("colour", "fill"),guide=guide_legend(keywidth=2,reverse=T))



########################################
# Descriptive stan model on each prompt
########################################

# Get the list of prompts:
Prompts = levels(Data$Prompt)

# Precompile the descriptive models (<2min):
binary_default <- stan_model(file = "StanModels/binary_default.stan",model_name="binary_default")
Likert_default <- stan_model(file = "StanModels/Likert_default.stan",model_name="Likert_default")


# Fit for each prompt and answer type (~5min for me, run fewer chains if you have less than 16 cores):
t <- Sys.time()
binary_fits <- list()
Likert_fits <- list()
for (prompt in Prompts){
  prompt_data <- Data %>%
    filter(Prompt==prompt & AnswerType=="Binary")
  stan_binary_model_data <- list(
    N = nrow(prompt_data),
    evidence_zero = if_else(prompt_data$Evidence==0,1,0),
    evidence_one = if_else(prompt_data$Evidence==1,1,0),
    scaled_evidence = prompt_data$scaled_evidence,
    truth = prompt_data$centered_truth,
    y = as.integer(prompt_data$NormAnswer)
  )
  prompt_data <- Data %>%
    filter(Prompt==prompt & AnswerType=="Likert")
  stan_Likert_model_data <- list(
    N = nrow(prompt_data),
    evidence_zero = if_else(prompt_data$Evidence==0,1,0),
    evidence_one = if_else(prompt_data$Evidence==1,1,0),
    scaled_evidence = prompt_data$scaled_evidence,
    truth = prompt_data$centered_truth,
    y = as.integer(prompt_data$Answer)
  )
  binary_fits[[prompt]] <- sampling(binary_default,
                                data = stan_binary_model_data,
                                chains = 14, iter = 4000,warmup=1000,
                                check_data=F,refresh=0
  )
  Likert_fits[[prompt]] <- sampling(Likert_default,
                                data = stan_Likert_model_data,
                                chains = 14, iter = 4000,warmup=1000,
                                # For Likert-True, the model requires a higher adapt_delta to avoid divergent transitions:
                                control=list(adapt_delta=ifelse(prompt=="True",0.95,0.8)),
                                check_data=F,refresh=0,include=F,pars="predictor"
  )
  if(sum(get_divergent_iterations(binary_fits[[prompt]]))>0){cat(sum(get_divergent_iterations(binary_fits[[prompt]]))," divergent transitions for ",prompt," binary model.\n")}
  if(sum(get_divergent_iterations(Likert_fits[[prompt]]))>0){cat(sum(get_divergent_iterations(Likert_fits[[prompt]]))," divergent transitions for ",prompt," Likert model.\n")}
  cat(paste0("Prompt ",prompt," finished.\n"))
}
print(Sys.time()-t)

# Generate a tibble with all parameters posterior estimates:
parameters_tibble = tibble(Prompt=character(),Response=character(),Parameter=character(),mean=numeric(),`2.5%`=numeric(),`97.5%`=numeric())
for (prompt in Prompts){
  tmp <- summary(binary_fits[[prompt]],pars = c("beta_truth","beta_ev","beta_inter","offset_z","offset_o"))$summary %>%
    as_tibble(rownames="Parameter") %>%
    mutate(Prompt = prompt,
           Response = "binary") %>%
    select(Prompt,Response,Parameter,mean,`2.5%`,`97.5%`)
  parameters_tibble <- bind_rows(parameters_tibble,tmp)
  tmp <- summary(Likert_fits[[prompt]],pars = c("beta_truth","beta_ev","beta_inter","offset_z","offset_o"))$summary %>%
    as_tibble(rownames="Parameter") %>%
    mutate(Prompt = prompt,
           Response = "Likert") %>%
    select(Prompt,Response,Parameter,mean,`2.5%`,`97.5%`)
  parameters_tibble <- bind_rows(parameters_tibble,tmp)
}

# Print a latex table for the paper (still requires some improvement):
parameters_tibble %>%
  xtable() %>%
  print(.,include.rownames=F)


#######################################
# Theory-driven models on each prompt
#######################################

# Precompile theory-driven models (~4min):
t <- Sys.time()
binary_truth <- stan_model(file = "StanModels/binary_truth.stan",model_name="binary_truth")
Likert_truth <- stan_model(file = "StanModels/Likert_truth.stan",model_name="Likert_truth")
binary_evidence <- stan_model(file = "StanModels/binary_evidence.stan",model_name="binary_evidence")
Likert_evidence <- stan_model(file = "StanModels/Likert_evidence.stan",model_name="Likert_evidence")
binary_assertability <- stan_model(file = "StanModels/binary_assertability.stan",model_name="binary_assertability")
Likert_assertability <- stan_model(file = "StanModels/Likert_assertability.stan",model_name="Likert_assertability")
Sys.time()-t

# Fit all models (~30min):
# We need a LOT of samples to reliably compute Bayes' factors
binary_truth_fits <- list()
Likert_truth_fits <- list()
binary_evidence_fits <- list()
Likert_evidence_fits <- list()
binary_assertability_fits <- list()
Likert_assertability_fits <- list()
t <- Sys.time()
print(t)
for (prompt in Prompts){
  prompt_data <- Data %>%
    filter(Prompt==prompt & AnswerType=="Binary")
  stan_binary_model_data <- list(
    N = nrow(prompt_data),
    evidence_zero = if_else(prompt_data$Evidence==0,1,0),
    evidence_one = if_else(prompt_data$Evidence==1,1,0),
    scaled_evidence = prompt_data$Evidence,
    truth = prompt_data$centered_truth,
    y = as.integer(prompt_data$NormAnswer)
  )
  prompt_data <- Data %>%
    filter(Prompt==prompt & AnswerType=="Likert")
  stan_Likert_model_data <- list(
    N = nrow(prompt_data),
    evidence_zero = if_else(prompt_data$Evidence==0,1,0),
    evidence_one = if_else(prompt_data$Evidence==1,1,0),
    scaled_evidence = prompt_data$Evidence,
    truth = prompt_data$centered_truth,
    y = as.integer(prompt_data$Answer)
  )
  # Lots of copy-pasting that could be simplified here:
  binary_truth_fits[[prompt]] <- sampling(binary_truth,
                                    data = stan_binary_model_data,
                                    chains = 14, iter = 20000,warmup=1000,
                                    check_data=F,refresh=0,control=list(adapt_delta=.9)
  )
  if(sum(get_divergent_iterations(binary_truth_fits[[prompt]]))>0){
    cat(sum(get_divergent_iterations(binary_truth_fits[[prompt]])),"divergent transitions for",prompt,"binary Truth model.\n")
  } else {
    cat(paste0(prompt," binary Truth model finished\n"))
  }
  binary_evidence_fits[[prompt]] <- sampling(binary_evidence,
                                          data = stan_binary_model_data,
                                          chains = 14, iter = 20000,warmup=1000,
                                          check_data=F,refresh=0,control=list(adapt_delta=.9)
  )
  if(sum(get_divergent_iterations(binary_evidence_fits[[prompt]]))>0){
    cat(sum(get_divergent_iterations(binary_evidence_fits[[prompt]])),"divergent transitions for",prompt,"binary Evidence model.\n")
  } else {
    cat(paste0(prompt," binary Evidence model finished\n"))
  }
  binary_assertability_fits[[prompt]] <- sampling(binary_assertability,
                                          data = stan_binary_model_data,
                                          chains = 14, iter = 20000,warmup=1000,
                                          check_data=F,refresh=0,control=list(adapt_delta=.9)
  )
  if(sum(get_divergent_iterations(binary_assertability_fits[[prompt]]))>0){
    cat(sum(get_divergent_iterations(binary_assertability_fits[[prompt]])),"divergent transitions for",prompt,"binary assertability model.\n")
  } else {
    cat(paste0(prompt," binary assertability model finished\n"))
  }
  Likert_truth_fits[[prompt]] <- sampling(Likert_truth,
                                          data = stan_Likert_model_data,
                                          chains = 14, iter = 20000,warmup=1000,
                                          check_data=F,refresh=0,
                                          # Increase adapt_delta for the True prompt to avoid divergent transition
                                          control=list(adapt_delta=ifelse(prompt=="True",.99,.9)),include=F,pars="predictor"
  )
  if(sum(get_divergent_iterations(Likert_truth_fits[[prompt]]))>0){
    cat(sum(get_divergent_iterations(Likert_truth_fits[[prompt]])),"divergent transitions for",prompt,"Likert Truth model.\n")
  } else {
    cat(paste0(prompt," Likert Truth model finished\n"))
  }
  Likert_evidence_fits[[prompt]] <- sampling(Likert_evidence,
                                             data = stan_Likert_model_data,
                                             chains = 14, iter = 20000,warmup=1000,
                                             check_data=F,refresh=0,control=list(adapt_delta=.9),include=F,pars="predictor"
  )
  if(sum(get_divergent_iterations(Likert_evidence_fits[[prompt]]))>0){
    cat(sum(get_divergent_iterations(Likert_evidence_fits[[prompt]])),"divergent transitions for",prompt,"Likert Evidence model.\n")
  } else {
    cat(paste0(prompt," Likert Evidence model finished\n"))
  }
  Likert_assertability_fits[[prompt]] <- sampling(Likert_assertability,
                                            data = stan_Likert_model_data,
                                            chains = 14, iter = 20000,warmup=1000,
                                            check_data=F,refresh=0,control=list(adapt_delta=.9),include=F,pars="predictor"
  )
  if(sum(get_divergent_iterations(Likert_assertability_fits[[prompt]]))>0){
    cat(sum(get_divergent_iterations(Likert_assertability_fits[[prompt]])),"divergent transitions for",prompt,"Likert assertability model.\n")
  } else {
    cat(paste0(prompt," Likert assertability model finished\n"))
  }
  cat(paste0("Prompt ",prompt," finished after ",round(Sys.time()-t),"min.\n"))
}


# Compute marginal log-likelihood (~15min):
marginal_loglik_binary <- list()
marginal_loglik_Likert <- list()
t <- Sys.time()
for (prompt in Prompts){
  tmp_binary_list <- list()
  tmp_Likert_list <- list()
  for (model in c("truth","evidence","assertability")){
    tmp_binary_list[[model]] <- eval(parse(text=paste0("bridge_sampler(binary_",model,"_fits[[\"",prompt,"\"]], silent = TRUE,cores = 16)")))
    tmp_Likert_list[[model]] <- eval(parse(text=paste0("bridge_sampler(Likert_",model,"_fits[[\"",prompt,"\"]], silent = TRUE,cores = 16)")))
  }
  marginal_loglik_binary[[prompt]] <- tmp_binary_list
  marginal_loglik_Likert[[prompt]] <- tmp_Likert_list
  print(paste0("Prompt ",prompt," finished."))
}
Sys.time()-t

# Function to generate clean report with Bayes factors:
report_BF <- function(MLL_list,model_type="Models"){
  tmp <- unlist(lapply(MLL_list,function(x)x$logml))
  or <- order(-tmp)
  bf12 <- as.numeric(bayes_factor(MLL_list[[or[1]]],MLL_list[[or[2]]])$bf)
  bf23 <- as.numeric(bayes_factor(MLL_list[[or[2]]],MLL_list[[or[3]]])$bf)
  cat(paste0(model_type," order: ",paste(c("truth","evidence","assertability")[or],collapse = " > "),"\n"))
  cat(paste0("Bayes factor in favor of first over second: ",signif(bf12,3),"\n"))
  cat(paste0("Bayes factor in favor of second over third: ",signif(bf23,3),"\n"))
}

# Apply to each prompt and answer type:
for (prompt in Prompts){
  cat(paste0("Prompt: ",prompt,"\n"))
  report_BF(marginal_loglik_binary[[prompt]],model_type="Binary models")
  report_BF(marginal_loglik_Likert[[prompt]],model_type="Likert models")
  cat("\n")
}



##############################################################
# Plot predictions of the Truth/Evidence/Assertability models
##############################################################


# Function to compute mean prediction and CI on the mean prediction from parameter samples for one combination of predictors:
# Note that it would be more efficient computation-wise to vectorize it along all predictors values,
# but it would probably use too much memory (~53M prediction points). Doing it rowwise works, and it just takes around 1min per model.
compute_binary_prediction <- function(scaled_evidence,truth,evidence_zero,evidence_one,params){
  Prediction <- with(params,
                     alpha
                     + beta_truth * truth
                     + beta_ev * scaled_evidence
                     + beta_inter * truth * scaled_evidence
                     - evidence_zero*offset_z
                     + evidence_one*offset_o)
  Prediction = plogis(Prediction)
  CI = ci(Prediction,method="HDI")
  return(tibble(mean_pred = mean(Prediction),
                ci_lo = CI$CI_low,
                ci_hi = CI$CI_high))
}

# For proportional odd logistic regression, one can show that the expected predicted value is:
# \sum XP(X) = 7 - P(X\leq 6)- ... - P(X\leq 1) = 7 - f(c6-pred) ...  - f(c1-pred) where f is the inv-logit function.

compute_Likert_prediction <- function(scaled_evidence,truth,evidence_zero,evidence_one,params){
  Predictor <- with(params,
                     beta_truth * truth
                     + beta_ev * scaled_evidence
                     + beta_inter * truth * scaled_evidence
                     - evidence_zero*offset_z
                     + evidence_one*offset_o)
  Prediction = 7 - rowSums(plogis(sweep(params$c,1,Predictor,"-")))
  CI = ci(Prediction,method="HDI")
  return(tibble(mean_pred = mean(Prediction),
                ci_lo = CI$CI_low,
                ci_hi = CI$CI_high))
}

cluster_copy(cluster,c("compute_binary_prediction","compute_Likert_prediction"))


# Sample of scaled evidence along which we'll generate predictions:
ScaledEvidence=seq(min(Data$Evidence),max(Data$Evidence),length.out=100L)

# Tibble with all combinations of predictors for which we'll generate predictions
predictions_template <- expand_grid(Prompt=levels(Data$Prompt),truth=c(-0.5,0.5),scaled_evidence=ScaledEvidence) %>%
  mutate(Evidence = (scaled_evidence-min(Data$Evidence))/(max(Data$Evidence-min(Data$Evidence))),
         evidence_zero = ifelse(Evidence==0,1,0),
         evidence_one = ifelse(Evidence==1,1,0)) %>%
  filter(!(truth==max(truth)&evidence_zero==1)&!(truth==min(truth)&evidence_one==1)) %>% # remove impossible combinations
  mutate(Prompt=factor(Prompt,levels=levels(Data$Prompt)))


# manual dodge between data and model prediction at min and max evidence:
dodge_amount=0.015

plotdata <- Data %>%
  mutate(Truth=if_else(Truth,"aTrue","aFalse"),
         Answer=as.numeric(Answer),
         Evidence = case_when(
           Evidence==0 ~ dodge_amount,
           Evidence==1 ~ 1+dodge_amount,
           T ~ Evidence
         )) %>%
  group_by(AnswerType,Prompt,Evidence,Truth) %>%
  summarize(mean=mean(Answer),
            se=se(Answer),
            y_min=mean-se,
            y_max=mean+se)

Likert_plotdata <- plotdata %>%
  filter(AnswerType=="Likert")

binary_plotdata <- plotdata %>%
  filter(AnswerType=="Binary")



# Functions to plot the result
plot_binary_predictions <- function(data){
  data %>%
    filter(Evidence>.05&Evidence<.95) %>%
    ggplot(aes(x=Evidence-dodge_amount*(evidence_zero+evidence_one),
               y=mean_pred,ymin=ci_lo,ymax=ci_hi,
               color=Truth,fill=Truth,shape=Truth))+
    facet_wrap(~Prompt)+
    #geom_line()+
    geom_ribbon(alpha=.2,color="transparent")+
    geom_point(data=filter(data,evidence_zero>0|evidence_one>0),size=2)+
    geom_errorbar(data=filter(data,evidence_zero>0|evidence_one>0),width=.05)+
    geom_point(data=binary_plotdata ,aes(x=Evidence,y=mean,col=Truth,shape=Truth),inherit.aes = F)+
    geom_errorbar(data=binary_plotdata ,aes(x=Evidence,y=mean,ymin=y_min,ymax=y_max,col=Truth),width=.05)+
    theme_bw()+
    ylab("Acceptability")+
    coord_cartesian(ylim=c(-0.01,1.01))+
    scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
    scale_color_manual(
      values=c("#F5793A","#0F2080","#E69F00", "#56B4E9"),
      name="Truth:",
      labels=c("False (data)","True (data)","False (model)","True (model)"),
      aesthetics = c("colour", "fill"),
      guide=guide_legend(keywidth=2,reverse=T))+
    scale_shape_manual(values=c(15,15,18,18),name="Truth:",labels=c("False (data)","True (data)","False (model)","True (model)"),guide=guide_legend(keywidth=2,reverse=T))
}
plot_Likert_predictions <- function(data){
  data %>%
    filter(Evidence>.05&Evidence<.95) %>%
    ggplot(aes(x=Evidence-dodge_amount*(evidence_zero+evidence_one),
               y=mean_pred,ymin=ci_lo,ymax=ci_hi,
               color=Truth,fill=Truth,shape=Truth))+
    facet_wrap(~Prompt)+
    geom_line()+
    geom_ribbon(alpha=.2,color="transparent")+
    geom_point(data=filter(data,evidence_zero>0|evidence_one>0),size=2)+
    geom_errorbar(data=filter(data,evidence_zero>0|evidence_one>0),width=.05)+
    geom_point(data=Likert_plotdata ,aes(x=Evidence,y=mean,col=Truth,shape=Truth),inherit.aes = F)+
    geom_errorbar(data=Likert_plotdata ,aes(x=Evidence,y=mean,ymin=y_min,ymax=y_max,col=Truth),width=.05)+
    theme_bw()+
    ylab("Acceptability")+
    coord_cartesian(ylim=c(.99,7.01))+
    scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
    scale_y_continuous(breaks=1:7,minor_breaks = NULL)+
    scale_color_manual(
      values=c("#F5793A","#0F2080","#E69F00", "#56B4E9"),
      name="Truth:",
      labels=c("False (data)","True (data)","False (model)","True (model)"),
      aesthetics = c("colour", "fill"),
      guide=guide_legend(keywidth=2,reverse=T))+
    scale_shape_manual(values=c(15,15,18,18),name="Truth:",labels=c("False (data)","True (data)","False (model)","True (model)"),guide=guide_legend(keywidth=2,reverse=T))
}



# Binary truth model

# Extract parameters for all prompts:
binary_truth_pars_samples <- lapply(as.list(Prompts),function(prompt)as_tibble(extract(binary_truth_fits[[prompt]],pars=c("alpha","beta_truth","beta_ev","beta_inter","offset_z","offset_o"))))
names(binary_truth_pars_samples) <- Prompts
cluster_copy(cluster,"binary_truth_pars_samples")

# Generate predictions
binary_truth_predictions <- predictions_template %>%
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    compute_binary_prediction(scaled_evidence,truth,evidence_zero,evidence_one,binary_truth_pars_samples[[Prompt]])
  ) %>%
  collect() %>%
  ungroup() %>%
  mutate(Truth=if_else(truth==max(truth),"True","False"))

# Plot
plot_binary_predictions(binary_truth_predictions)

# Binary evidence model

# Extract parameters for all prompts:
binary_evidence_pars_samples <- lapply(as.list(Prompts),function(prompt)as_tibble(extract(binary_evidence_fits[[prompt]],pars=c("alpha","beta_truth","beta_ev","beta_inter","offset_z","offset_o"))))
names(binary_evidence_pars_samples) <- Prompts
cluster_copy(cluster,"binary_evidence_pars_samples")

# Generate predictions
binary_evidence_predictions <- predictions_template %>%
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    compute_binary_prediction(scaled_evidence,truth,evidence_zero,evidence_one,binary_evidence_pars_samples[[Prompt]])
  ) %>%
  collect() %>%
  ungroup() %>%
  mutate(Truth=if_else(truth==max(truth),"True","False"))

# Plot
plot_binary_predictions(binary_evidence_predictions)

# Binary assertability model

# Extract parameters for all prompts:
binary_assertability_pars_samples <- lapply(as.list(Prompts),function(prompt)as_tibble(extract(binary_assertability_fits[[prompt]],pars=c("alpha","beta_truth","beta_ev","beta_inter","offset_z","offset_o"))))
names(binary_assertability_pars_samples) <- Prompts
cluster_copy(cluster,"binary_assertability_pars_samples")

# Generate predictions
binary_assertability_predictions <- predictions_template %>%
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    compute_binary_prediction(scaled_evidence,truth,evidence_zero,evidence_one,binary_assertability_pars_samples[[Prompt]])
  ) %>%
  collect() %>%
  ungroup() %>%
  mutate(Truth=if_else(truth==max(truth),"True","False"))

# Plot
plot_binary_predictions(binary_assertability_predictions)


# Likert truth model

# Extract parameters for all prompts:
Likert_truth_pars_samples <- lapply(as.list(Prompts),function(prompt)as_tibble(extract(Likert_truth_fits[[prompt]],pars=c("c","beta_truth","beta_ev","beta_inter","offset_z","offset_o"))))
names(Likert_truth_pars_samples) <- Prompts
cluster_copy(cluster,"Likert_truth_pars_samples")

# Generate predictions
Likert_truth_predictions <- predictions_template %>%
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    compute_Likert_prediction(scaled_evidence,truth,evidence_zero,evidence_one,Likert_truth_pars_samples[[Prompt]])
  ) %>%
  collect() %>%
  ungroup() %>%
  mutate(Truth=if_else(truth==max(truth),"True","False"))

# Plot
plot_Likert_predictions(Likert_truth_predictions)

# Likert evidence model

# Extract parameters for all prompts:
Likert_evidence_pars_samples <- lapply(as.list(Prompts),function(prompt)as_tibble(extract(Likert_evidence_fits[[prompt]],pars=c("c","beta_truth","beta_ev","beta_inter","offset_z","offset_o"))))
names(Likert_evidence_pars_samples) <- Prompts
cluster_copy(cluster,"Likert_evidence_pars_samples")

# Generate predictions
Likert_evidence_predictions <- predictions_template %>%
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    compute_Likert_prediction(scaled_evidence,truth,evidence_zero,evidence_one,Likert_evidence_pars_samples[[Prompt]])
  ) %>%
  collect() %>%
  ungroup() %>%
  mutate(Truth=if_else(truth==max(truth),"True","False"))

# Plot
plot_Likert_predictions(Likert_evidence_predictions)


# Likert assertability model

# Extract parameters for all prompts:
Likert_assertability_pars_samples <- lapply(as.list(Prompts),function(prompt)as_tibble(extract(Likert_assertability_fits[[prompt]],pars=c("c","beta_truth","beta_ev","beta_inter","offset_z","offset_o"))))
names(Likert_assertability_pars_samples) <- Prompts
cluster_copy(cluster,"Likert_assertability_pars_samples")

# Generate predictions
Likert_assertability_predictions <- predictions_template %>%
  rowwise() %>%
  partition(cluster) %>%
  mutate(
    compute_Likert_prediction(scaled_evidence,truth,evidence_zero,evidence_one,Likert_assertability_pars_samples[[Prompt]])
  ) %>%
  collect() %>%
  ungroup() %>%
  mutate(Truth=if_else(truth==max(truth),"True","False"))

# Plot
plot_Likert_predictions(Likert_assertability_predictions)



##################################
# Plot best model for each prompt
##################################

find_best_model <- function(prompt,answer_type){
  if(answer_type=="Likert"){
    MLL_list <- marginal_loglik_Likert
  } else {
      MLL_list <- marginal_loglik_binary
  }
  lapply(MLL_list[[prompt]],function(x)x$logml) |>
  unlist() |> which.max() |> names()
}

best_model_tibble <- expand_grid(Prompt=factor(Prompts,levels=Prompts),AnswerType=c("Binary","Likert")) %>%
  rowwise() %>%
  mutate(BestModel = find_best_model(Prompt,AnswerType)) %>%
  ungroup()


binary_best_predictions <- binary_truth_predictions %>%
  inner_join(binary_evidence_predictions,
             by=c("Prompt","Truth","truth","scaled_evidence","Evidence","evidence_zero","evidence_one"),
             suffix=c("_truth","_evidence")) %>%
  inner_join(binary_assertability_predictions,
             by=c("Prompt","Truth","truth","scaled_evidence","Evidence","evidence_zero","evidence_one"),
             suffix=c("","_assertability")) %>%
  rename_with(~paste0(.,"_assertability"),c(mean_pred,ci_lo,ci_hi)) %>%
  left_join(filter(best_model_tibble,AnswerType=="Binary")) %>%
  mutate(
    mean_pred = case_when(
      BestModel=="truth" ~ mean_pred_truth,
      BestModel=="evidence" ~ mean_pred_evidence,
      BestModel=="assertability" ~ mean_pred_assertability
    ),
    ci_lo = case_when(
      BestModel=="truth" ~ ci_lo_truth,
      BestModel=="evidence" ~ ci_lo_evidence,
      BestModel=="assertability" ~ ci_lo_assertability
    ),
    ci_hi = case_when(
      BestModel=="truth" ~ ci_hi_truth,
      BestModel=="evidence" ~ ci_hi_evidence,
      BestModel=="assertability" ~ ci_hi_assertability
    )
  ) %>%
  select(-starts_with("mean_pred_"),-starts_with("ci_lo_"),-starts_with("ci_hi_")) %>%
  mutate(Facet=paste0(Prompt,": ",BestModel," model"))

tmp <- filter(best_model_tibble,AnswerType=="Binary") %>%
  mutate(Facet=paste0(Prompt,": ",BestModel," model"))
binary_facet_labels <- tmp$Facet
names(binary_facet_labels) = Prompts

# Figure 3a:
p <- binary_best_predictions %>%
  filter(Evidence>.05&Evidence<.95) %>%
  ggplot(aes(x=Evidence-dodge_amount*(evidence_zero+evidence_one),
             y=mean_pred,ymin=ci_lo,ymax=ci_hi,
             color=Truth,fill=Truth,shape=Truth))+
  facet_wrap(~Prompt,labeller = labeller(Prompt=binary_facet_labels))+
  geom_line()+
  geom_ribbon(alpha=.2,color="transparent")+
  geom_point(data=filter(binary_best_predictions,evidence_zero>0|evidence_one>0),size=2)+
  geom_errorbar(data=filter(binary_best_predictions,evidence_zero>0|evidence_one>0),width=.05)+
  geom_point(data=binary_plotdata ,aes(x=Evidence,y=mean,col=Truth,shape=Truth),inherit.aes = F)+
  geom_errorbar(data=binary_plotdata ,aes(x=Evidence,y=mean,ymin=y_min,ymax=y_max,col=Truth),width=.05)+
  theme_bw()+
  ylab("Acceptability")+
  coord_cartesian(ylim=c(-0.01,1.01))+
  scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
  scale_color_manual(
    values=c("#F5793A","#0F2080","#E69F00", "#56B4E9"),
    name="Truth:",
    labels=c("False (data)","True (data)","False (model)","True (model)"),
    aesthetics = c("colour", "fill"),
    guide=guide_legend(keywidth=2,reverse=T))+
  scale_shape_manual(values=c(15,15,18,18),name="Truth:",labels=c("False (data)","True (data)","False (model)","True (model)"),guide=guide_legend(keywidth=2,reverse=T))

pdf(file="Graphs/Best_Predictions_binary.pdf",width=9,height=5)
print(p)
dev.off()


Likert_best_predictions <- Likert_truth_predictions %>%
  inner_join(Likert_evidence_predictions,
             by=c("Prompt","Truth","truth","scaled_evidence","Evidence","evidence_zero","evidence_one"),
             suffix=c("_truth","_evidence")) %>%
  inner_join(Likert_assertability_predictions,
             by=c("Prompt","Truth","truth","scaled_evidence","Evidence","evidence_zero","evidence_one"),
             suffix=c("","_assertability")) %>%
  rename_with(~paste0(.,"_assertability"),c(mean_pred,ci_lo,ci_hi)) %>%
  left_join(filter(best_model_tibble,AnswerType=="Likert")) %>%
  mutate(
    mean_pred = case_when(
      BestModel=="truth" ~ mean_pred_truth,
      BestModel=="evidence" ~ mean_pred_evidence,
      BestModel=="assertability" ~ mean_pred_assertability
    ),
    ci_lo = case_when(
      BestModel=="truth" ~ ci_lo_truth,
      BestModel=="evidence" ~ ci_lo_evidence,
      BestModel=="assertability" ~ ci_lo_assertability
    ),
    ci_hi = case_when(
      BestModel=="truth" ~ ci_hi_truth,
      BestModel=="evidence" ~ ci_hi_evidence,
      BestModel=="assertability" ~ ci_hi_assertability
    )
  ) %>%
  select(-starts_with("mean_pred_"),-starts_with("ci_lo_"),-starts_with("ci_hi_")) %>%
  mutate(Facet=paste0(Prompt,": ",BestModel," model"))

tmp <- filter(best_model_tibble,AnswerType=="Likert") %>%
  mutate(Facet=paste0(Prompt,": ",BestModel," model"))
Likert_facet_labels <- tmp$Facet
names(Likert_facet_labels) = Prompts

# Figure 3b:
p <- Likert_best_predictions %>%
  filter(Evidence>.05&Evidence<.95) %>%
  ggplot(aes(x=Evidence-dodge_amount*(evidence_zero+evidence_one),
             y=mean_pred,ymin=ci_lo,ymax=ci_hi,
             color=Truth,fill=Truth,shape=Truth))+
  facet_wrap(~Prompt,labeller = labeller(Prompt=Likert_facet_labels))+
  geom_line()+
  geom_ribbon(alpha=.2,color="transparent")+
  geom_point(data=filter(Likert_best_predictions,evidence_zero>0|evidence_one>0),size=2)+
  geom_errorbar(data=filter(Likert_best_predictions,evidence_zero>0|evidence_one>0),width=.05)+
  geom_point(data=Likert_plotdata ,aes(x=Evidence,y=mean,col=Truth,shape=Truth),inherit.aes = F)+
  geom_errorbar(data=Likert_plotdata ,aes(x=Evidence,y=mean,ymin=y_min,ymax=y_max,col=Truth),width=.05)+
  theme_bw()+
  ylab("Acceptability")+
  coord_cartesian(ylim=c(.99,7.01))+
  scale_y_continuous(breaks=1:7,minor_breaks = NULL)+
  scale_x_continuous(name="Evidence",breaks = seq(0,1,.25),labels=round(seq(0,1,.25),2))+
  scale_color_manual(
    values=c("#F5793A","#0F2080","#E69F00", "#56B4E9"),
    name="Truth:",
    labels=c("False (data)","True (data)","False (model)","True (model)"),
    aesthetics = c("colour", "fill"),
    guide=guide_legend(keywidth=2,reverse=T))+
  scale_shape_manual(values=c(15,15,18,18),name="Truth:",labels=c("False (data)","True (data)","False (model)","True (model)"),guide=guide_legend(keywidth=2,reverse=T))

pdf(file="Graphs/Best_Predictions_Likert.pdf",width=9,height=5)
print(p)
dev.off()


#######################
# Effect size measures
#######################

# We use the ordinal superiority measure discussed in Ryu&Agresti (2008) as our effect size.
# gamma = P(Y1>Y2) + 0.5P(Y1=Y2)
# Most importantly, this measure is valid for both Likert and binary.
# We can also extract credible intervals on these effect sizes.


# General function for arbitrary Likert scales:
Likert_ordinal_superiority = function(eta1,eta2,c_vec,k){
  # c_vec must be length k-1
  c_vec=c(-Inf,c_vec,Inf)
  sum((plogis(eta2-c_vec[1:k])-plogis(eta2-c_vec[2:(k+1)]))*plogis(eta1-c_vec[2:(k+1)])+
        0.5*(plogis(eta2-c_vec[1:k])-plogis(eta2-c_vec[2:(k+1)]))*(plogis(eta1-c_vec[1:k])-plogis(eta1-c_vec[2:(k+1)])))
}
# More efficient here:
Likert_ordinal_superiority = function(eta1,eta2,c1,c2,c3,c4,c5,c6){
  # This is stupid, but not doing it generates useless warning which slows things down A LOT:
  eta1=c(eta1)
  eta2=c(eta2)
  c_vec=c(-Inf,c1,c2,c3,c4,c5,c6,Inf)
  sum((plogis(eta2-c_vec[1:7])-plogis(eta2-c_vec[2:8]))*plogis(eta1-c_vec[2:8])+
        0.5*(plogis(eta2-c_vec[1:7])-plogis(eta2-c_vec[2:8]))*(plogis(eta1-c_vec[1:7])-plogis(eta1-c_vec[2:8])))
}
binary_ordinal_superiority = function(eta1,eta2,beta0){
  0.5*(1+plogis(eta1+beta0)-plogis(eta2+beta0))
}
min_evid=min(Data$scaled_evidence)
max_evid=max(Data$scaled_evidence)

cluster_copy(cluster,c("binary_ordinal_superiority","Likert_ordinal_superiority","min_evid","max_evid"))

binary_effect_size <- function(model){
  tmp_parameters <- extract(model,pars = c("beta_truth","beta_ev","beta_inter","offset_z","offset_o","alpha"))
  parameters <- bind_cols(
    tibble(beta_truth=tmp_parameters[[1]]),
    tibble(beta_ev=tmp_parameters[[2]]),
    tibble(beta_inter=tmp_parameters[[3]]),
    tibble(offset_z=tmp_parameters[[4]]),
    tibble(offset_o=tmp_parameters[[5]]),
    tibble(alpha=tmp_parameters[[6]])
  ) %>%
    rowwise() %>%
    partition(cluster) %>%
    mutate(
      ES_truth = binary_ordinal_superiority(eta1=0.5*beta_truth,eta2=-0.5*beta_truth,beta0=alpha),
      ES_ev = binary_ordinal_superiority(eta1=0.5*beta_ev,eta2=-0.5*beta_ev,beta0=alpha),
      ES_inter = binary_ordinal_superiority(eta1=0.25*beta_inter,eta2=-0.25*beta_inter,beta0=alpha),
      ES_zero = binary_ordinal_superiority(eta1=-0.5*beta_truth+min_evid*beta_ev-0.5*min_evid*beta_inter+offset_z,
                                           eta2=-0.5*beta_truth+min_evid*beta_ev-0.5*min_evid*beta_inter,beta0=alpha),
      ES_one = binary_ordinal_superiority(eta1=0.5*beta_truth+max_evid*beta_ev+0.5*max_evid*beta_inter+offset_o,
                                          eta2=0.5*beta_truth+max_evid*beta_ev+0.5*max_evid*beta_inter,beta0=alpha)
    ) %>%
    collect() %>%
    ungroup()
  parameters %>%
    summarise(across(starts_with("ES_"),list(mean=mean,CIlo=~ci(.x,method="HDI")$CI_low,CIhi=~ci(.x,method="HDI")$CI_high))) %>%
    pivot_longer(starts_with("ES_"),names_to = c(NA,"parameter",".value"),names_sep = "_")
}
Likert_effect_size <- function(model){
  tmp_parameters <- extract(model,pars = c("beta_truth","beta_ev","beta_inter","offset_z","offset_o","c"))
  parameters <- bind_cols(
    tibble(beta_truth=tmp_parameters[[1]]),
    tibble(beta_ev=tmp_parameters[[2]]),
    tibble(beta_inter=tmp_parameters[[3]]),
    tibble(offset_z=tmp_parameters[[4]]),
    tibble(offset_o=tmp_parameters[[5]]),
    as_tibble(tmp_parameters[[6]]) %>% set_names(paste0("c",1:6))
  ) %>%
    rowwise() %>%
    partition(cluster) %>%
    mutate(
      ES_truth = Likert_ordinal_superiority(eta1=0.5*beta_truth,eta2=-0.5*beta_truth,c1,c2,c3,c4,c5,c6),
      ES_ev = Likert_ordinal_superiority(eta1=0.5*beta_ev,eta2=-0.5*beta_ev,c1,c2,c3,c4,c5,c6),
      ES_inter = Likert_ordinal_superiority(eta1=0.25*beta_inter,eta2=-0.25*beta_inter,c1,c2,c3,c4,c5,c6),
      ES_zero = Likert_ordinal_superiority(eta1=-0.5*beta_truth+min_evid*beta_ev-0.5*min_evid*beta_inter+offset_z,
                                           eta2=-0.5*beta_truth+min_evid*beta_ev-0.5*min_evid*beta_inter,c1,c2,c3,c4,c5,c6),
      ES_one = Likert_ordinal_superiority(eta1=0.5*beta_truth+max_evid*beta_ev+0.5*max_evid*beta_inter+offset_o,
                                           eta2=0.5*beta_truth+max_evid*beta_ev+0.5*max_evid*beta_inter,c1,c2,c3,c4,c5,c6)
    ) %>%
    collect() %>%
    ungroup()
  parameters %>%
    summarise(across(starts_with("ES_"),list(mean=mean,CIlo=~ci(.x,method="HDI")$CI_low,CIhi=~ci(.x,method="HDI")$CI_high))) %>%
    pivot_longer(starts_with("ES_"),names_to = c(NA,"parameter",".value"),names_sep = "_")
}

binary_effect_sizes <- list()
Likert_effect_sizes <- list()

# ~ 3min
t <- Sys.time()
for(prompt in Prompts){
  binary_effect_sizes[[prompt]] <- binary_effect_size(binary_fits[[prompt]])
  Likert_effect_sizes[[prompt]] <- Likert_effect_size(Likert_fits[[prompt]])
  cat("Prompt ",prompt," done in ",Sys.time()-t,".\n")
}

# saveRDS(binary_fits,"binary_fits.RDS")
# saveRDS(Likert_fits,"Likert_fits.RDS")

binary_fits <- readRDS("binary_fits.RDS")
Likert_fits <- readRDS("Likert_fits.RDS")


effect_sizes_tibble = tibble(Prompt=character(),Response=character(),Parameter=character(),mean=numeric(),`2.5%`=numeric(),`97.5%`=numeric())
for (prompt in Prompts){
  tmp <- binary_effect_sizes[[prompt]] %>%
    mutate(Prompt = prompt,
           Response = "binary") %>%
    select(Prompt,Response,everything()) %>%
    set_names(c("Prompt","Response","Parameter","mean","2.5%","97.5%"))
  effect_sizes_tibble <- bind_rows(effect_sizes_tibble,tmp)
  tmp <- Likert_effect_sizes[[prompt]] %>%
    mutate(Prompt = prompt,
           Response = "Likert") %>%
    select(Prompt,Response,everything()) %>%
    set_names(c("Prompt","Response","Parameter","mean","2.5%","97.5%"))
  effect_sizes_tibble <- bind_rows(effect_sizes_tibble,tmp)
}

effect_sizes_tibble <- effect_sizes_tibble %>%
  rename(gamma=mean) %>%
  mutate(Parameter=case_when(
    Parameter=="truth" ~ "truth",
    Parameter=="ev" ~ "evidence",
    Parameter=="inter" ~ "interaction",
    Parameter=="zero" ~ "offset_z",
    Parameter=="one" ~ "offset_o"
  )) %>%
  mutate(
    gamma_sym = if_else(gamma>0.5,gamma,1-gamma),
    min_sample_size = ifelse(Response=="binary",loess_min_sample_binary(gamma_sym),loess_min_sample_Likert(gamma_sym)),
    min_sample_size = pmax(3,ceiling(min_sample_size))
  ) %>%
  select(-gamma_sym)

# Add gamma to parameters table for comparison:
parameters_tibble %>%
  left_join(select(effect_sizes_tibble,Prompt,Response,Parameter,gamma)) %>%
  mutate(gamma=round(gamma,2)) %>%
  print(n=Inf)

# Print a latex table for the paper:
effect_sizes_tibble %>%
  xtable() %>%
  print(.,include.rownames=F)

