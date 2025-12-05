rm(list=ls())

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/simulation/param_realdata.zeisel.imputed.RData")

keep <- with(
  as.data.frame(match_params),
  kon  > 1  & kon  < 10 &
    koff > 1 & koff < 10 &
    s    > 20  & s    < 1000
)
match_params <- match_params[keep, , drop = FALSE]

# Optional: slightly boost s
match_params[, "s"] <- match_params[, "s"] * 2

save(match_params,
     file = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/yenchi-jerry_yuhong/out/veloUncertainty/kevin/simulation/param_realdata.zeisel.imputed_KZL.RData")
