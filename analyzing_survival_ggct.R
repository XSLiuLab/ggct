load("./ggct_survival_data_on_metastasis_and_stages.Rdata")
library(survival)
library(survminer)

# cnv_sm and exp_sm is for metastasis (pathological_M)
# cnv_sm2 and expr_sm2 is for stages(pathological_stage (TNM))

metastasis <- ggct_on_metastasis[, c(-4,-5)]
stages <- ggct_on_stages[, c(-4, -5)]

groupSruvival <- function(df, event="OS_IND", time="OS", var=NULL, time.limit=NULL, interval=c("open","close"),
                          method=c("quartile", "mean", "median", "percent", "custom"), percent=NULL, 
                          step=20, custom_fun=NULL, group1="High", group2="Low"){
    #'@param df a data.frame at leaset including three column which refer to survival info and variable used to
    #'                       set group.
    #'@param event The status indicator, normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). 
    #'For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
    #' Although unusual, the event indicator can be omitted, in which case all subjects are assumed to have an event.
    #'@param time for right censored data, this is the follow up time. For interval data, the first argument is the starting time for the interval.
    #'@param var a column name of df data.frame which specify the numeric column used to group. 
    #'@param time.limit a number, if any time bigger than it, the corresponding row will be ruled out from df data.frame
    #'@param interval specify if the border should be included. 
    #'@param method a mehod used to group, you can use one of "quartile", "mean", "median", "percent", "custom". 
    #'@param percent if you choose "pecent" method, you must specify what percent you want.
    #'@param step a number to define a step length for finding p values and time cutoff.
    #'@param custom_fun a function object. If you choose "custom" method, you must define a function which use
    #'values of 'var' column as input and return a list with two logical vectors responding group1 and group2, respectively.                       
    #'@param group1 a character define the group name which respond to bigger values in quartile, mean, median
    #',percent methods and the first logical vector in custom method. 'High' is default.
    #'@param group2 a character define the group name which respond to smaller values in quartile, mean, median
    #',percent methods and the second logical vector in custom method. 'Low' is default.
    
    if(is.null(var)){
        stop("Please set your variable name for building groups!")
    }    
    
    if(! "data.frame" %in% class(df)){
        stop("Be care about you input data. The df argument should be a data.frame!")
    }
    
    method=match.arg(method)
    
    if(method == "percent" & is.null(percent)){
        stop("You choose 'percent' method but not percent values defined in argument. Please check!")
    }
    
    df1 <- df[, c(event, time, var)]
    
    if(!is.null(time.limit)){
        df1 <- df[df[,time]<=time.limit,]
    }
    
    if (method == "custom"){
        if(is.null(custom_fun)){
            stop("You should specify a function to group 'var' column.")
        }
        if(!inherits(custom_fun, "function")){
            stop("custom_fun should be a function type.")
        }
        th <- custom_fun(as.numeric(df1[, var]))
        if(typeof(th) != "list" | length(th)!=2){
            stop("You should store your results for grouping in a list with two logical elements!")
        }
        if(any(th[[1]] & th[[2]])){
            stop("You are mapping one record to two different groups, please check your function!")
        }
        df1$group <- NA
        df1$group[th[[1]]] <- group1 
        df1$group[th[[2]]] <- group2
        df2 <- subset(df1, group %in% c(group1, group2))
        #df2$group[df2[,var]>th2] = "High"
    }else{
    
        sm <- summary(df1[, var])
        
        # choose threshold by method
        if (method == "quartile"){
            th1 <- as.numeric(sm[2])
            th2 <- as.numeric(sm[5])
        }
        if (method == "mean"){
            th1 <- th2 <- as.numeric(sm[4])
        }
        if (method == "median"){
            th1 <- th2 <- as.numeric(sm[3])
        }
        if (method == "percent"){
            th1 <- as.numeric(quantile(as.numeric(df1[, var]) ,percent))
            th2 <- as.numeric(quantile(as.numeric(df1[, var]) ,1-percent))
        }
    
        # do groupping
        interval <- match.arg(interval)
        
    
        if (interval == "open"){
            df2 <- df1[df[, var]>th2 | df[, var] < th1 , ]
            df2$group <- NA
            df2$group[df2[,var]>th2] <- group1
            df2$group[df2[,var]<th1] <- group2
        }else{
            df2 <- df1[df[, var]>=th2 | df[, var] <= th1 , ]
            df2$group <- NA
            df2$group[df2[,var]>=th2] <- group1
            df2$group[df2[,var]<=th1] <- group2
        }
    }
    
    colnames(df2)[colnames(df2) == time] <- "time"
    colnames(df2)[colnames(df2) == event] <- "event"
    sfit <- survfit(Surv(time, event) ~ group, data=df2)
    
    res <- list()
    res$data <- df2 
    
    search_cutoff <- function(mat, step=20){
        # mat <- os_mat
        days <- max(na.omit(mat$time))
        
        cutoff_list <- {}
        pval_list <- {}
        cutoff <- list()
        
        while(days > 300){
            if(length(table(mat$group))!=2){
                return(cutoff)
            }
            mat <- mat[mat$time<=days,]
            sdf <- survdiff(Surv(time, event)~group,data = mat)
            p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
            cutoff_list <- c(cutoff_list, days)
            pval_list <- c(pval_list, p.val)
            days <- days - step
        }
        
        min_pval <- min(pval_list)
        min_cutoff <- cutoff_list[which(pval_list==min_pval)]
        
        cutoff$cutoff_list <- cutoff_list
        cutoff$pval_list <- pval_list
        cutoff$min_pval <- min_pval
        cutoff$min_cutoff <- min_cutoff
        
        return(cutoff)
        
    }
    
    cutoff <- search_cutoff(df2, step=step)
    
    res$cutoff <- cutoff
    
    return(res)
} 


plot_surv <- function(os_mat, cutoff=NULL, pval=TRUE, ...){
    if(is.null(cutoff)){
        cutoff <- max(na.omit(os_mat$time))
    }
    fit <- survfit(Surv(time, event) ~ group, data = os_mat[os_mat$time<=cutoff,])
    ggsurv <- ggsurvplot(
        fit,                     # survfit object with calculated statistics.
        data = os_mat[os_mat$time<=cutoff,],           # data used to fit survival curves.
        risk.table = FALSE,       # show risk table.
        pval = pval,             # show p-value of log-rank test.
        conf.int = FALSE,         # show confidence intervals for
        # point estimates of survival curves.
        palette = c("red", "blue"),
        # palette = c("#E7B800", "#2E9FDF"),
        # xlim = c(0,5000),         # present narrower X axis, but not affect
        # survival estimates.
        xlab = "",   # customize X axis label.
        ylab = "",
        # break.time.by = 100,     # break X axis in time intervals by 500.
        # ggtheme = theme_light(), # customize plot and risk table with a theme.
        risk.table.y.text.col = T,# colour risk table text annotations.
        risk.table.height = 0.25, # the height of the risk table
        risk.table.y.text = FALSE,# show bars instead of names in text annotations
        # in legend of risk table.
        ncensor.plot = FALSE,      # plot the number of censored subjects at time t
        ncensor.plot.height = 0.25,
        # conf.int.style = "step",  # customize style of confidence intervals
        # surv.median.line = "hv",  # add the median survival pointer.
        legend.labs =
            c("High expression", "Low expression"),    # change legend labels.
        legend.title = "",
        size = 0.5,
        censor.size = 2,
        ...
    )
    
    ggsurv <- ggpar(
        ggsurv,
        font.caption = c(6, "plain", "black"),
        font.tickslab = c(8, "plain", "black"),
        font.legend = c(6, "plain", "black"),
        ggtheme = theme_survminer()+theme(legend.background=element_rect(fill=NA), legend.key.height = unit(3,"mm"),
                                          line = element_line(size=.1), legend.position =c(0.8, 0.9))
    )
    
    ggsurv$plot
}


defineGroups <- function(x){
    res <- list()
    res[[1]] <- x > 0.4
    res[[2]] <- x > -0.1 & x < 0.1
    return(res)
}

## test
test <- groupSruvival(df=ggct_on_stages, event="OS_IND", time="OS", 
                      var="cnv",time.limit = NULL, interval = "open", method = "custom",
                      custom_fun = defineGroups,
                      percent = NULL, 
                      step=5)

## do group survival analysis

## cnv effect on survival of metastasis
res.non_meta.cnv <- groupSruvival(df=subset(metastasis, pathologic_M=="M0"), event="OS_IND", time="OS", 
                                  var="cnv",time.limit = NULL, interval = "open", method = "custom",
                                  custom_fun = defineGroups, step=5)
head(res.non_meta.cnv$cutoff$pval_list[order(res.non_meta.cnv$cutoff$pval_list)],100)
head(res.non_meta.cnv$cutoff$cutoff_list[order(res.non_meta.cnv$cutoff$pval_list)],100)
plot_surv(os_mat = res.non_meta.cnv$data)

res.meta.cnv <- groupSruvival(df=subset(metastasis, pathologic_M=="M1"), event="OS_IND", time="OS", 
                              var="cnv",time.limit = NULL, interval = "open", method = "custom",
                              custom_fun = defineGroups, step=5)

head(res.meta.cnv$cutoff$pval_list[order(res.meta.cnv$cutoff$pval_list)],10)
head(res.meta.cnv$cutoff$cutoff_list[order(res.meta.cnv$cutoff$pval_list)],10)
plot_surv(os_mat = res.meta.cnv$data, cutoff = 3650)


## expression effect on survival of metastasis
res.non_meta.exp <- groupSruvival(df=subset(metastasis, pathologic_M=="M0"), event="OS_IND", time="OS", 
                                  var="expression",time.limit = NULL, interval = "open", method = "quartile", 
                                  percent = NULL, 
                                  step=5)
head(res.non_meta.exp$cutoff$pval_list[order(res.non_meta.exp$cutoff$pval_list)],100)
head(res.non_meta.exp$cutoff$cutoff_list[order(res.non_meta.exp$cutoff$pval_list)],100)
plot_surv(os_mat = res.non_meta.exp$data)

res.meta.exp <- groupSruvival(df=subset(metastasis, pathologic_M=="M1"), event="OS_IND", time="OS", 
                              var="expression",time.limit = NULL, interval = "open", method = "quartile", 
                              percent = NULL, 
                              step=5)

head(res.meta.exp$cutoff$pval_list[order(res.meta.exp$cutoff$pval_list)],10)
head(res.meta.exp$cutoff$cutoff_list[order(res.meta.exp$cutoff$pval_list)],10)
plot_surv(os_mat = res.meta.exp$data)




test <- groupSruvival(df=ggct_on_metastasis, event="OS_IND", time="OS", var="expression",time.limit = NULL,
                      interval = "open", method = "quartile", percent = NULL, pval=TRUE, pval.method=TRUE,
                      step=5)
plot_surv(os_mat = test$data, cutoff = test$cutoff$min_cutoff)






##
# write.table(test$data[test$data$time<=822,], file="./test_survival.txt", col.names = T, row.names = F, quote=F,
#             sep="\t")

# sfit <- survfit(Surv(df2[,"OS"], df2[, "OS_IND"])~group, data=df2)
# ggsurvplot(sfit)




# load("./df2.Rdata")
# test_fun <- function(df2, time="OS", event="OS_IND"){
#     sfit <- survfit(as.formula("Surv(OS, OS_IND)~group"), data=df2)
#     ggsurvplot(sfit, data=df2)
# }
# df1 <- df2; rm(df2)
# test_fun(df2=df1)
# 
# 
# 
# plot(test)
# ggsurvplot(test, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
#            legend.labs=c("High expression", "low"), legend.title="Metastasis status",  
#            palette=c("dodgerblue2", "orchid2"), 
#            main="Kaplan-Meier Curve for Lung Cancer Survival", 
#            risk.table.height=.15)
