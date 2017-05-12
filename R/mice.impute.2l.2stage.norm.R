mice.impute.2l.2stage.norm <-
function(y, ry, x,type,method_est="mm",...)
{
  return(mice.impute.2l.2stage.norm.intern(y, ry, x,type,method_est=method_est))
}
