#
#
#
#
#
#

# 增加 drop=F 的參數, 只取一欄的時候, 還能保留 data.frame 結構
my_get_variable=function (physeq, varName)
{
  if (is.null(sample_data(physeq, FALSE))) {
    stop("Your phyloseq data object does not have a sample-data component\n",
         "Try ?sample_data for more details.")
  }
  return(as(sample_data(physeq), "data.frame")[, varName,drop=F])
}
