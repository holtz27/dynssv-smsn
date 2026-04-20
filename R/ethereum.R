# Obtenha os dados do ethereum
ethereum = quantmod::getSymbols('ETH-USD', 
                               src='yahoo', 
                               from='2017-07-01', to='2022-12-06',
                               #from='2017-01-07', to='2022-10-12',
                               auto.assign=FALSE)
ethereum = na.omit(ethereum)
ethereum = data.frame(ethereum)
dates = as.Date(row.names(ethereum), '%Y-%m-%d')
ethereum = ethereum[,'ETH.USD.Adjusted']
#View(ethereum)
T = length(ethereum)
log.ret = 100*(log( ethereum[2:T])-log(ethereum[1:(T-1)]))
T = length(log.ret)
log.ret=log.ret-mean(log.ret)
#ytrain=log.ret

# Plots
library(ggplot2)
df = data.frame(Return=log.ret, Tempo=dates[-1])
g = ggplot(df) + geom_line(aes(x=Tempo, y=Return))
g = g + scale_x_date(date_breaks="10 month", date_labels="%b %Y")
g = g + theme_test() + theme(axis.title.y=element_text(size=18),
                             axis.text.x=element_text(size=11),
                             axis.text.y=element_text(size=18))
g = g + xlab('')
h = ggplot(df, aes(Return))
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color='white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x=element_text(size=18),
                             axis.text.x=element_text(size=12),
                             axis.text.y=element_text(size=18))
gridExtra::grid.arrange(g, h, nrow=1, ncol=2) 
data_summary = matrix(c( T, mean(log.ret),
                         sd(log.ret),
                         min(log.ret),
                         max(log.ret),
                         moments::skewness(log.ret),
                         moments::kurtosis(log.ret)), nrow=1)
colnames(data_summary)=c('T', 'mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round(data_summary, digits=3)

# Test
# Obtenha os dados do ethereum
ethereum = quantmod::getSymbols('ETH-USD', 
                               src='yahoo', 
                               #from='2017-01-07', to='2022-06-12',
                               from='2022-06-12', to='2022-10-12',
                               auto.assign=FALSE)
ethereum = na.omit(ethereum)
ethereum = data.frame(ethereum)
dates = as.Date(row.names(ethereum), '%Y-%m-%d')
ethereum = ethereum[,'ETH.USD.Adjusted']
#View(ethereum)
T = length(ethereum)
log.ret = 100*(log( ethereum[2:T])-log(ethereum[1:(T-1)]))
T = length(log.ret)
#log.ret=log.ret-mean(log.ret)
ytest=log.ret

y=c(ytrain, ytest)
