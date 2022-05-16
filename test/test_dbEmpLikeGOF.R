library('dbEmpLikeGOF')
x <- rnorm(10, 0,1)
dbEmpLikeGOF(x = x, testcall = "normal", pvl.Table = FALSE, num.mc = 5000)
