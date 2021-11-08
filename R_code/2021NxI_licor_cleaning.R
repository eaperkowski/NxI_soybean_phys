library(readLicorData)
library(data.table)

## Manually clean licor a/ci files
alb.oct14 <- licorData(location = "../licor_data/2021-10-14_aci_soybean_alb.xlsx")
alb.oct15 <- licorData(location = "../licor_data/2021-10-15_aci_soybean_alb.xlsx")
nvg.oct14 <- licorData(location = "../licor_data/2021-10-14_aci_soybean_nvg.xlsx")
nvg.oct15 <- licorData(location = "../licor_data/2021-10-15_aci_soybean_nvg.xlsx")
ozz.oct14 <- licorData(location = "../licor_data/2021-10-14_aci_soybean_ozz.xlsx")
ozz.oct15 <- licorData(location = "../licor_data/2021-10-15_aci_soybean_ozz.xlsx")
gib.oct15 <- licorData(location = "../licor_data/2021-10-15_aci_soybean_gib.xlsx")

## Save licor a/ci files
write.csv(alb.oct14, "../licor_data_cleaned/aci/2021-10-14_aci_soybean_alb.csv")
write.csv(alb.oct15, "../licor_data_cleaned/aci/2021-10-15_aci_soybean_alb.csv")
write.csv(nvg.oct14, "../licor_data_cleaned/aci/2021-10-14_aci_soybean_nvg.csv")
write.csv(nvg.oct15, "../licor_data_cleaned/aci/2021-10-15_aci_soybean_nvg.csv")
write.csv(ozz.oct14, "../licor_data_cleaned/aci/2021-10-14_aci_soybean_ozz.csv")
write.csv(ozz.oct15, "../licor_data_cleaned/aci/2021-10-15_aci_soybean_ozz.csv")
write.csv(gib.oct15, "../licor_data_cleaned/aci/2021-10-15_aci_soybean_gib.csv")

## Manually clean licor resp files
alb.oct14 <- licorData(location = "../licor_data/2021-10-14_resp_soybean_alb.xlsx")
alb.oct15 <- licorData(location = "../licor_data/2021-10-15_resp_soybean_alb.xlsx")
nvg.oct14 <- licorData(location = "../licor_data/2021-10-14_resp_soybean_nvg.xlsx")
nvg.oct15 <- licorData(location = "../licor_data/2021-10-15_resp_soybean_nvg.xlsx")
ozz.oct14 <- licorData(location = "../licor_data/2021-10-14_resp_soybean_ozz.xlsx")
ozz.oct15 <- licorData(location = "../licor_data/2021-10-15_resp_soybean_ozz.xlsx")
gib.oct15 <- licorData(location = "../licor_data/2021-10-15_resp_soybean_gib.xlsx")


## Save licor resp files
write.csv(alb.oct14, "../licor_data_cleaned/resp/2021-10-14_resp_soybean_alb.csv")
write.csv(alb.oct15, "../licor_data_cleaned/resp/2021-10-15_resp_soybean_alb.csv")
write.csv(nvg.oct14, "../licor_data_cleaned/resp/2021-10-14_resp_soybean_nvg.csv")
write.csv(nvg.oct15, "../licor_data_cleaned/resp/2021-10-15_resp_soybean_nvg.csv")
write.csv(ozz.oct14, "../licor_data_cleaned/resp/2021-10-14_resp_soybean_ozz.csv")
write.csv(ozz.oct15, "../licor_data_cleaned/resp/2021-10-15_resp_soybean_ozz.csv")
write.csv(gib.oct15, "../licor_data_cleaned/resp/2021-10-15_resp_soybean_gib.csv")
