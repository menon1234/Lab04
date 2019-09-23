linreg <- setRefClass("linreg",fields = list(formula = "formula",data = "data.frame"), methods = list(
linreg = function(){

}, print = function(){
  return("hi")
}
))

linreg$methods(plot = function(){
  return("this")
})
myObj <- linreg$new()
myObj$plot()
