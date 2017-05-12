mice.impute.2l.2stage.pmm <-
function(y, ry, x,type,method_est="mm",incluster=FALSE,k=5,...)
{
  return(mice.impute.2l.2stage.norm.intern(y, ry, x,type,method_est=method_est,method_draw="pmm",incluster=incluster,kk=k))
}
