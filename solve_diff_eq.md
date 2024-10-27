# Intro to Differential Equations
H Elliott

- [Basic model of population growth for a single
  species](#basic-model-of-population-growth-for-a-single-species)
- [If r = 0](#if-r-0)

## Basic model of population growth for a single species

Let ![N(t)](https://latex.codecogs.com/svg.latex?N%28t%29 "N(t)") be a
function that describes the size of a population evolving over time. Let
parameter ![r](https://latex.codecogs.com/svg.latex?r "r") describe the
“intrinsic” per-capita growth rate of the population - i.e., the
per-capita growth rate if the population was free from any resource
constraints or any other effects of its own density. Let parameter
![s](https://latex.codecogs.com/svg.latex?s "s") determine how the
current size of the population effects the per-capita growth rate.

We can write this as a first-order differential equation:

![\frac{1}{N} \frac{dN}{dt} = r + sN](https://latex.codecogs.com/svg.latex?%5Cfrac%7B1%7D%7BN%7D%20%5Cfrac%7BdN%7D%7Bdt%7D%20%3D%20r%20%2B%20sN "\frac{1}{N} \frac{dN}{dt} = r + sN")

We can solve this equation by separating variables:

- ![\frac{dN}{N ( r + s N)} = dt](https://latex.codecogs.com/svg.latex?%5Cfrac%7BdN%7D%7BN%20%28%20r%20%2B%20s%20N%29%7D%20%3D%20dt "\frac{dN}{N ( r + s N)} = dt")  
- ![\int \frac{dN}{N ( r + s N)} = \int dt](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7BdN%7D%7BN%20%28%20r%20%2B%20s%20N%29%7D%20%3D%20%5Cint%20dt "\int \frac{dN}{N ( r + s N)} = \int dt")
  - ![\int \frac{dN}{N ( r + s N)} = \int \left( \frac{1}{r} \frac{dN}{N} - \frac{1}{r} \frac{dN}{r + sN} \right)](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7BdN%7D%7BN%20%28%20r%20%2B%20s%20N%29%7D%20%3D%20%5Cint%20%5Cleft%28%20%5Cfrac%7B1%7D%7Br%7D%20%5Cfrac%7BdN%7D%7BN%7D%20-%20%5Cfrac%7B1%7D%7Br%7D%20%5Cfrac%7BdN%7D%7Br%20%2B%20sN%7D%20%5Cright%29 "\int \frac{dN}{N ( r + s N)} = \int \left( \frac{1}{r} \frac{dN}{N} - \frac{1}{r} \frac{dN}{r + sN} \right)")
    (using partial fractions)
  - ![= \frac{1}{r} \int \frac{dN}{N} - \frac{1}{r} \int \frac{dN}{r + sN}](https://latex.codecogs.com/svg.latex?%3D%20%5Cfrac%7B1%7D%7Br%7D%20%5Cint%20%5Cfrac%7BdN%7D%7BN%7D%20-%20%5Cfrac%7B1%7D%7Br%7D%20%5Cint%20%5Cfrac%7BdN%7D%7Br%20%2B%20sN%7D "= \frac{1}{r} \int \frac{dN}{N} - \frac{1}{r} \int \frac{dN}{r + sN}")  
  - ![= \frac{1}{r} ( \ln N - \ln (r + sN) )](https://latex.codecogs.com/svg.latex?%3D%20%5Cfrac%7B1%7D%7Br%7D%20%28%20%5Cln%20N%20-%20%5Cln%20%28r%20%2B%20sN%29%20%29 "= \frac{1}{r} ( \ln N - \ln (r + sN) )")  
- ![\frac{1}{r} (\ln N - \ln (Ns + r)) = t + C](https://latex.codecogs.com/svg.latex?%5Cfrac%7B1%7D%7Br%7D%20%28%5Cln%20N%20-%20%5Cln%20%28Ns%20%2B%20r%29%29%20%3D%20t%20%2B%20C "\frac{1}{r} (\ln N - \ln (Ns + r)) = t + C")  
- ![\ln N - \ln (Ns + r) = rt + Cr](https://latex.codecogs.com/svg.latex?%5Cln%20N%20-%20%5Cln%20%28Ns%20%2B%20r%29%20%3D%20rt%20%2B%20Cr "\ln N - \ln (Ns + r) = rt + Cr")  
- ![\ln \frac{N}{Ns + r} = rt + Cr](https://latex.codecogs.com/svg.latex?%5Cln%20%5Cfrac%7BN%7D%7BNs%20%2B%20r%7D%20%3D%20rt%20%2B%20Cr "\ln \frac{N}{Ns + r} = rt + Cr")  
- ![\frac{N}{Ns + r} = e^{rt + Cr}](https://latex.codecogs.com/svg.latex?%5Cfrac%7BN%7D%7BNs%20%2B%20r%7D%20%3D%20e%5E%7Brt%20%2B%20Cr%7D "\frac{N}{Ns + r} = e^{rt + Cr}")  
- ![N = e^{rt + Cr} (Ns + r)](https://latex.codecogs.com/svg.latex?N%20%3D%20e%5E%7Brt%20%2B%20Cr%7D%20%28Ns%20%2B%20r%29 "N = e^{rt + Cr} (Ns + r)")  
- ![N = e^{rt + Cr} Ns + e^{rt + Cr} r](https://latex.codecogs.com/svg.latex?N%20%3D%20e%5E%7Brt%20%2B%20Cr%7D%20Ns%20%2B%20e%5E%7Brt%20%2B%20Cr%7D%20r "N = e^{rt + Cr} Ns + e^{rt + Cr} r")  
- ![N - e^{rt + Cr} Ns = e^{rt + Cr} r](https://latex.codecogs.com/svg.latex?N%20-%20e%5E%7Brt%20%2B%20Cr%7D%20Ns%20%3D%20e%5E%7Brt%20%2B%20Cr%7D%20r "N - e^{rt + Cr} Ns = e^{rt + Cr} r")  
- ![N (1 - e^{rt + Cr} s) = e^{rt + Cr} r](https://latex.codecogs.com/svg.latex?N%20%281%20-%20e%5E%7Brt%20%2B%20Cr%7D%20s%29%20%3D%20e%5E%7Brt%20%2B%20Cr%7D%20r "N (1 - e^{rt + Cr} s) = e^{rt + Cr} r")  
- ![N = \frac{e^{rt + Cr} r}{1 - e^{rt + Cr} s}](https://latex.codecogs.com/svg.latex?N%20%3D%20%5Cfrac%7Be%5E%7Brt%20%2B%20Cr%7D%20r%7D%7B1%20-%20e%5E%7Brt%20%2B%20Cr%7D%20s%7D "N = \frac{e^{rt + Cr} r}{1 - e^{rt + Cr} s}")  
- ![N = \frac{r}{e^{-rt - Cr} - s}](https://latex.codecogs.com/svg.latex?N%20%3D%20%5Cfrac%7Br%7D%7Be%5E%7B-rt%20-%20Cr%7D%20-%20s%7D "N = \frac{r}{e^{-rt - Cr} - s}")  
- ![N = \frac{r}{A e^{-rt} - s}](https://latex.codecogs.com/svg.latex?N%20%3D%20%5Cfrac%7Br%7D%7BA%20e%5E%7B-rt%7D%20-%20s%7D "N = \frac{r}{A e^{-rt} - s}"),
  where
  ![A = e^{-Cr}](https://latex.codecogs.com/svg.latex?A%20%3D%20e%5E%7B-Cr%7D "A = e^{-Cr}")
  is a constant.

If at
![t = 0](https://latex.codecogs.com/svg.latex?t%20%3D%200 "t = 0"),
![N = N_0](https://latex.codecogs.com/svg.latex?N%20%3D%20N_0 "N = N_0")
(some initial population size), then:

- ![N_0 = \frac{r}{A e^{-r \cdot 0} - s}](https://latex.codecogs.com/svg.latex?N_0%20%3D%20%5Cfrac%7Br%7D%7BA%20e%5E%7B-r%20%5Ccdot%200%7D%20-%20s%7D "N_0 = \frac{r}{A e^{-r \cdot 0} - s}")  
- ![N_0 = \frac{r}{A - s}](https://latex.codecogs.com/svg.latex?N_0%20%3D%20%5Cfrac%7Br%7D%7BA%20-%20s%7D "N_0 = \frac{r}{A - s}")  
- ![N_0 A - N_0 s = r](https://latex.codecogs.com/svg.latex?N_0%20A%20-%20N_0%20s%20%3D%20r "N_0 A - N_0 s = r")  
- ![A = \frac{r + N_0 s}{N_0}](https://latex.codecogs.com/svg.latex?A%20%3D%20%5Cfrac%7Br%20%2B%20N_0%20s%7D%7BN_0%7D "A = \frac{r + N_0 s}{N_0}")

So the particular solution, for any initial population
![N_0](https://latex.codecogs.com/svg.latex?N_0 "N_0") at time
![t = 0](https://latex.codecogs.com/svg.latex?t%20%3D%200 "t = 0"):

![N(t) = \frac{r}{\frac{r + N_0 s}{N_0} e^{-rt} - s}](https://latex.codecogs.com/svg.latex?N%28t%29%20%3D%20%5Cfrac%7Br%7D%7B%5Cfrac%7Br%20%2B%20N_0%20s%7D%7BN_0%7D%20e%5E%7B-rt%7D%20-%20s%7D "N(t) = \frac{r}{\frac{r + N_0 s}{N_0} e^{-rt} - s}")

Note that solving the equation relied on
![r \neq 0](https://latex.codecogs.com/svg.latex?r%20%5Cneq%200 "r \neq 0")
since we saw ![1/r](https://latex.codecogs.com/svg.latex?1%2Fr "1/r")
appear early in the process. Further, examining the final solution you
can see that
![r = 0](https://latex.codecogs.com/svg.latex?r%20%3D%200 "r = 0")
implies that
![N = 0 / (s - s) = 0/0](https://latex.codecogs.com/svg.latex?N%20%3D%200%20%2F%20%28s%20-%20s%29%20%3D%200%2F0 "N = 0 / (s - s) = 0/0")
which is undefined.

## If r = 0

Notice that if
![r = 0](https://latex.codecogs.com/svg.latex?r%20%3D%200 "r = 0"), then
the differential equation simplifies so that population growth is
entirely dependent on the current population size - there is no
intrinsic growth:

![\frac{dN}{dt} = sN^2](https://latex.codecogs.com/svg.latex?%5Cfrac%7BdN%7D%7Bdt%7D%20%3D%20sN%5E2 "\frac{dN}{dt} = sN^2")

- ![\frac{dN}{N^2} = s dt](https://latex.codecogs.com/svg.latex?%5Cfrac%7BdN%7D%7BN%5E2%7D%20%3D%20s%20dt "\frac{dN}{N^2} = s dt")  
- ![\int \frac{1}{N^2} dN = \int s dt](https://latex.codecogs.com/svg.latex?%5Cint%20%5Cfrac%7B1%7D%7BN%5E2%7D%20dN%20%3D%20%5Cint%20s%20dt "\int \frac{1}{N^2} dN = \int s dt")  
- ![- \frac{1}{N} = st + C](https://latex.codecogs.com/svg.latex?-%20%5Cfrac%7B1%7D%7BN%7D%20%3D%20st%20%2B%20C "- \frac{1}{N} = st + C")  
- ![N = \frac{1}{-st - C}](https://latex.codecogs.com/svg.latex?N%20%3D%20%5Cfrac%7B1%7D%7B-st%20-%20C%7D "N = \frac{1}{-st - C}")

If
![N = N_0](https://latex.codecogs.com/svg.latex?N%20%3D%20N_0 "N = N_0")
at ![t = 0](https://latex.codecogs.com/svg.latex?t%20%3D%200 "t = 0"),
then:

- ![N_0 = \frac{1}{-s \cdot 0 - C}](https://latex.codecogs.com/svg.latex?N_0%20%3D%20%5Cfrac%7B1%7D%7B-s%20%5Ccdot%200%20-%20C%7D "N_0 = \frac{1}{-s \cdot 0 - C}")  
- ![N_0 = \frac{1}{-C}](https://latex.codecogs.com/svg.latex?N_0%20%3D%20%5Cfrac%7B1%7D%7B-C%7D "N_0 = \frac{1}{-C}")  
- ![C = -\frac{1}{N_0}](https://latex.codecogs.com/svg.latex?C%20%3D%20-%5Cfrac%7B1%7D%7BN_0%7D "C = -\frac{1}{N_0}")

So the particular solution is:

![N(t) = \frac{1}{-s t + \frac{1}{N_0}}](https://latex.codecogs.com/svg.latex?N%28t%29%20%3D%20%5Cfrac%7B1%7D%7B-s%20t%20%2B%20%5Cfrac%7B1%7D%7BN_0%7D%7D "N(t) = \frac{1}{-s t + \frac{1}{N_0}}")

Now it’s clear that this solution will fail if
![N_0 = 0](https://latex.codecogs.com/svg.latex?N_0%20%3D%200 "N_0 = 0"),
but that is already a special case, since if
![N_0 = 0](https://latex.codecogs.com/svg.latex?N_0%20%3D%200 "N_0 = 0")
and ![r = 0](https://latex.codecogs.com/svg.latex?r%20%3D%200 "r = 0"),
the the population will remain at
![0](https://latex.codecogs.com/svg.latex?0 "0") for all time. In this
case,
![N(t) = 0](https://latex.codecogs.com/svg.latex?N%28t%29%20%3D%200 "N(t) = 0")
would be the solution.
