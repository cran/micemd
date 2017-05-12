mice.impute.2l.2stage.bin <-
function(y, ry, x,type,method_est="mm",...)
{
  return(mice.impute.2l.2stage.bin.intern(y=y, ry=ry, x=x,type=type,method_est=method_est))
}
