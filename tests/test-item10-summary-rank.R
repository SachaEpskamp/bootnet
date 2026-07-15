# Regression test for audit item 10: fix the rank-handling typo cluster in
# summary.bootnet() (R/summaryMethod.R).
#
# The sampleTable value_min column was never populated: the rank = TRUE branch
# self-assigned rank_min and the rank = FALSE branch wrote to a stray rank_min
# column instead of value_min. bootTable and sampleTable are now handled
# symmetrically.
#
# Downstream summary computations consume only bootTable's value_min/value_max,
# so the DEFAULT summary output (rank = FALSE) must be IDENTICAL pre/post fix.
# The reference below was captured from the pre-fix HEAD build (same seed/data).

library(bootnet)

set.seed(12345)
data <- matrix(rnorm(200 * 5), ncol = 5)
net <- estimateNetwork(data, default = "pcor")
boot <- bootnet(net, nBoots = 5, nCores = 1, type = "nonparametric")

## --- Default summary (rank = FALSE): must match pre-fix output exactly --------
s1 <- summary(boot)
stopifnot(is.data.frame(s1))

s1df <- as.data.frame(s1)
num <- round(s1df[, sapply(s1df, is.numeric), drop = FALSE], 8)

ref <- structure(list(sample = c(0.03449199, 0.05519022, 0.0227629,
-0.02307723, 0.06255573, -0.00520549, 0.02004666, -0.17685221,
0.07231285, 0.0708044, 0.13552235, 0.12229988, 0.36691101, 0.27562499,
0.18624114), mean = c(0.09034526, 0.01905111, -0.04415811, 0.04986248,
0.03803585, -0.0128202, -0.00592234, -0.15316589, 0.07536826,
0.13462012, 0.25266972, 0.24716407, 0.29747807, 0.43837887, 0.30959486
), sd = c(0.0393745, 0.03680817, 0.09622281, 0.05856796, 0.0687994,
0.10007649, 0.07525739, 0.1018841, 0.03718448, 0.04474448, 0.09414054,
0.14950656, 0.13090777, 0.18500741, 0.1560004), CIlower = c(-0.04425701,
-0.01842611, -0.16968273, -0.14021315, -0.07504307, -0.20535848,
-0.13046812, -0.3806204, -0.0020561, -0.01868456, -0.05275874,
-0.17671325, 0.10509547, -0.09438983, -0.12575966), CIupper = c(0.113241,
0.12880656, 0.21520853, 0.0940587, 0.20015453, 0.1949475, 0.17056145,
0.02691599, 0.14668181, 0.16029336, 0.32380343, 0.421313, 0.62872656,
0.64563981, 0.49824195), q2.5 = c(0.03600311, -0.01693724, -0.14220502,
-0.0072907, -0.00362341, -0.16994109, -0.11677137, -0.3156466,
0.02330554, 0.08919065, 0.17765781, 0.09762702, 0.12524752, 0.19263502,
0.18308938), q97.5 = c(0.1323914, 0.05996744, 0.08349734, 0.12498373,
0.15895771, 0.09716566, 0.0947911, -0.03286105, 0.12779181, 0.20314961,
0.3880692, 0.48330586, 0.48924441, 0.68530716, 0.57269652), q2.5_non0 = c(0.03600311,
-0.01693724, -0.14220502, -0.0072907, -0.00362341, -0.16994109,
-0.11677137, -0.3156466, 0.02330554, 0.08919065, 0.17765781,
0.09762702, 0.12524752, 0.19263502, 0.18308938), mean_non0 = c(0.09034526,
0.01905111, -0.04415811, 0.04986248, 0.03803585, -0.0128202,
-0.00592234, -0.15316589, 0.07536826, 0.13462012, 0.25266972,
0.24716407, 0.29747807, 0.43837887, 0.30959486), q97.5_non0 = c(0.1323914,
0.05996744, 0.08349734, 0.12498373, 0.15895771, 0.09716566, 0.0947911,
-0.03286105, 0.12779181, 0.20314961, 0.3880692, 0.48330586, 0.48924441,
0.68530716, 0.57269652), var_non0 = c(0.00155035, 0.00135484,
0.00925883, 0.00343021, 0.00473336, 0.0100153, 0.00566368, 0.01038037,
0.00138269, 0.00200207, 0.00886244, 0.02235221, 0.01713684, 0.03422774,
0.02433613), sd_non0 = c(0.0393745, 0.03680817, 0.09622281, 0.05856796,
0.0687994, 0.10007649, 0.07525739, 0.1018841, 0.03718448, 0.04474448,
0.09414054, 0.14950656, 0.13090777, 0.18500741, 0.1560004), prop0 = c(0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), row.names = c(NA,
-15L), class = "data.frame")

stopifnot(identical(names(num), names(ref)))
stopifnot(nrow(num) == nrow(ref))
for (nm in names(ref)) {
  stopifnot(isTRUE(all.equal(num[[nm]], ref[[nm]], tolerance = 1e-8)))
}

## --- rank = TRUE completes and returns a data frame ---------------------------
s2 <- summary(boot, rank = TRUE)
stopifnot(is.data.frame(s2))

cat("test-item10-summary-rank.R: OK\n")
