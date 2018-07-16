# CAPM-credit

Simple Markovitz (mean-variance) optimization on a bond portfolio, plus an application of random matrix theory on the correlation matrix estimation.

I wrote a detailed report on this work which you can find on my website antoine-falck.fr.

## Data

The data used come from Bloomberg, they gather bond prices of the index `Markit iBoxxEUR High Yield Corporates BB Top 50 Mid Price TCA TRI`, from 2016-01-01 to 2018-03-29. 

## Mean-variance optimization

We want to maximize the yield's expectation of our portfolio under a variance constraint, based on the past observations. We discuss the case with and without a non-risky asset.

You can refer to the report I wrote or to *H. Markowitz,Portfolio Selection, The Journal of Finance, 7 (1952), pp. 77–91* for the complete reference.

## Random matrix theory

The previous section implies the estimation of the correlation matrix. Based on the random matrix theory we can make a new "filtered" estimation of this matrix, which gives us new efficient frontiers.

Again you can refer to the report or fin the complete reference in *L. Laloux, P. Cizeau, J.-P. Bouchaud, and M. Potters,Noise Dressing ofFinancial Correlation Matrices, Physical Review Letters, 83 (1999), pp. 1467–1470*.

