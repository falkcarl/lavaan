### Carl F. Falk
### Last update: 2018-05-09
### Test of "Robust" difference tests

context("lav_test_LRT")

###################################################################
## Test satorra.bentler.2001, satorra.bentler.2010 and satorra.2000
## Without starting values
###################################################################

HS.model <- '
    visual  =~ x1 + lam2*x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
'

## Default settings

## unrestricted model
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
              auto.cov.lv.x=TRUE, estimator="MLM",
              meanstructure=TRUE, int.ov.free=TRUE)

## restricted model
const2<-"lam2 == -1.03"
fit.const <- lavaan(HS.model, data=HolzingerSwineford1939,constraints=const2,
                    auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
                    auto.cov.lv.x=TRUE, estimator="MLM",meanstructure=TRUE, int.ov.free=TRUE)

lavSB2001<-lavTestLRT(fit,fit.const,method="satorra.bentler.2001")
lavSB2010<-lavTestLRT(fit,fit.const,method="satorra.bentler.2010")
lavS2000<-lavTestLRT(fit,fit.const,method="satorra.2000",A.method="exact")

test_that("Testing robust difference tests, default (lavaan)", {
  expect_equal(198.37, lavSB2001$`Chisq diff`[2], tolerance=.1)
  expect_equal(302.09, lavSB2010$`Chisq diff`[2], tolerance=.1)
  expect_equal(261.62, lavS2000$`Chisq diff`[2], tolerance=.1)
})

## mimic="Mplus"

## unrestricted model
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
              auto.cov.lv.x=TRUE, estimator="MLM",
              meanstructure=TRUE, int.ov.free=TRUE, mimic="Mplus")

## restricted model
const2<-"lam2 == -1.03"
fit.const <- lavaan(HS.model, data=HolzingerSwineford1939,constraints=const2,
                    auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
                    auto.cov.lv.x=TRUE, estimator="MLM",meanstructure=TRUE, int.ov.free=TRUE,
                    mimic="Mplus")

mplusSB2001<-lavTestLRT(fit,fit.const,method="satorra.bentler.2001")
mplusSB2010<-lavTestLRT(fit,fit.const,method="satorra.bentler.2010")
mplusS2000<-lavTestLRT(fit,fit.const,method="satorra.2000",A.method="exact")

test_that("Testing robust difference tests, Mplus", {
  expect_equal(171.81, mplusSB2001$`Chisq diff`[2], tolerance=.1)
  expect_equal(256.69, mplusSB2010$`Chisq diff`[2], tolerance=.1)
  expect_equal(261.62, mplusS2000$`Chisq diff`[2], tolerance=.1)
})

## mimic="EQS"

## unrestricted model
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
              auto.cov.lv.x=TRUE, estimator="MLM",
              meanstructure=TRUE, int.ov.free=TRUE, mimic="EQS")

## restricted model
const2<-"lam2 == -1.03"
fit.const <- lavaan(HS.model, data=HolzingerSwineford1939,constraints=const2,
                    auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
                    auto.cov.lv.x=TRUE, estimator="MLM",meanstructure=TRUE, int.ov.free=TRUE,
                    mimic="EQS")

eqsSB2001<-lavTestLRT(fit,fit.const,method="satorra.bentler.2001")
eqsSB2010<-lavTestLRT(fit,fit.const,method="satorra.bentler.2010")
eqsS2000<-lavTestLRT(fit,fit.const,method="satorra.2000",A.method="exact")

test_that("Testing robust difference tests, EQS", {
  expect_equal(198.54, eqsSB2001$`Chisq diff`[2], tolerance=.1)
  expect_equal(302.46, eqsSB2010$`Chisq diff`[2], tolerance=.1)
  expect_equal(262.12, eqsS2000$`Chisq diff`[2], tolerance=.1)
})

###################################################################
## Test satorra.bentler.2001, satorra.bentler.2010 and satorra.2000
## With some starting values, but tests should be same as above.
## sb2010 in particular can break if wrong starting values are used
## for model M10
###################################################################

HS.model <- ' visual =~ x1 + start(0.7)*lam2*x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9
x1 ~~ start(0.5)*x1 + psi1*x1
'

## Default settings

## unrestricted model
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
              auto.cov.lv.x=TRUE, estimator="MLM",
              meanstructure=TRUE, int.ov.free=TRUE)

## restricted model
const2<-"lam2 == -1.03"
fit.const <- lavaan(HS.model, data=HolzingerSwineford1939,constraints=const2,
                    auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
                    auto.cov.lv.x=TRUE, estimator="MLM",meanstructure=TRUE, int.ov.free=TRUE)

lavSB2001<-lavTestLRT(fit,fit.const,method="satorra.bentler.2001")
lavSB2010<-lavTestLRT(fit,fit.const,method="satorra.bentler.2010")
lavS2000<-lavTestLRT(fit,fit.const,method="satorra.2000",A.method="exact")

test_that("Testing robust difference tests, default (lavaan)", {
  expect_equal(198.37, lavSB2001$`Chisq diff`[2], tolerance=.1)
  expect_equal(302.09, lavSB2010$`Chisq diff`[2], tolerance=.1)
  expect_equal(261.62, lavS2000$`Chisq diff`[2], tolerance=.1)
})

## mimic="Mplus"

## unrestricted model
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
              auto.cov.lv.x=TRUE, estimator="MLM",
              meanstructure=TRUE, int.ov.free=TRUE, mimic="Mplus")

## restricted model
const2<-"lam2 == -1.03"
fit.const <- lavaan(HS.model, data=HolzingerSwineford1939,constraints=const2,
                    auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
                    auto.cov.lv.x=TRUE, estimator="MLM",meanstructure=TRUE, int.ov.free=TRUE,
                    mimic="Mplus")

mplusSB2001<-lavTestLRT(fit,fit.const,method="satorra.bentler.2001")
mplusSB2010<-lavTestLRT(fit,fit.const,method="satorra.bentler.2010")
mplusS2000<-lavTestLRT(fit,fit.const,method="satorra.2000",A.method="exact")

test_that("Testing robust difference tests, Mplus", {
  expect_equal(171.81, mplusSB2001$`Chisq diff`[2], tolerance=.1)
  expect_equal(256.69, mplusSB2010$`Chisq diff`[2], tolerance=.1)
  expect_equal(261.62, mplusS2000$`Chisq diff`[2], tolerance=.1)
})

## mimic="EQS"

## unrestricted model
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
              auto.cov.lv.x=TRUE, estimator="MLM",
              meanstructure=TRUE, int.ov.free=TRUE, mimic="EQS")

## restricted model
const2<-"lam2 == -1.03"
fit.const <- lavaan(HS.model, data=HolzingerSwineford1939,constraints=const2,
                    auto.var=TRUE, auto.fix.first=FALSE, std.lv=TRUE,
                    auto.cov.lv.x=TRUE, estimator="MLM",meanstructure=TRUE, int.ov.free=TRUE,
                    mimic="EQS")

eqsSB2001<-lavTestLRT(fit,fit.const,method="satorra.bentler.2001")
eqsSB2010<-lavTestLRT(fit,fit.const,method="satorra.bentler.2010")
eqsS2000<-lavTestLRT(fit,fit.const,method="satorra.2000",A.method="exact")

test_that("Testing robust difference tests, EQS", {
  expect_equal(198.54, eqsSB2001$`Chisq diff`[2], tolerance=.1)
  expect_equal(302.46, eqsSB2010$`Chisq diff`[2], tolerance=.1)
  expect_equal(262.12, eqsS2000$`Chisq diff`[2], tolerance=.1)
})

