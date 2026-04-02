if (!("future.apply" %in% rownames(installed.packages()))) {
  install.packages("future.apply")
}
library(future.apply)

if (!("parallel" %in% rownames(installed.packages()))) {
  install.packages("parallel")
}
library(parallel)

if (!("whitening" %in% rownames(installed.packages()))) {
  install.packages("whitening")
}
library(whitening)

if (!("lmerTest" %in% rownames(installed.packages()))) {
  install.packages("lmerTest")
}
library(lmerTest)

if (!("Rmpfr" %in% rownames(installed.packages()))) {
  install.packages("Rmpfr")
}
library(Rmpfr)

if (!("Matrix" %in% rownames(installed.packages()))) {
  install.packages("Matrix")
}
library(Matrix)

if (!("performance" %in% rownames(installed.packages()))) {
  install.packages("performance")
}
library(performance)

if (!("mvtnorm" %in% rownames(installed.packages()))) {
  install.packages("mvtnorm")
}
library(mvtnorm)

if (!("CR2" %in% rownames(installed.packages()))) {
  install.packages("CR2")
}
library(CR2)

if (!("clubSandwich" %in% rownames(installed.packages()))) {
  install.packages("clubSandwich")
}
library(clubSandwich)

if (!("ggplot2" %in% rownames(installed.packages()))) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!("dplyr" %in% rownames(installed.packages()))) {
  install.packages("dplyr")
}
library(dplyr)

if (!("purrr" %in% rownames(installed.packages()))) {
  install.packages("purrr")
}
library(purrr)

if (!("stringr" %in% rownames(installed.packages()))) {
  install.packages("stringr")
}
library(stringr)
