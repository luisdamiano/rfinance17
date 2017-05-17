# A Quick Intro to Hidden Markov Models Applied to Stock Volatility

Both the [presentation](./presentation/) and the [notebook](./notebook/) are part of the material presented in [R/Finance](http://www.rinfinance.com/) 2017.

## About R/Finance
Applied Finance with R From the inaugural conference in 2009, the annual R/Finance conference in Chicago has become the primary meeting for academics and practioners interested in using R in Finance. Participants from academia and industry mingle for two days to exchange ideas about current research, best practices and applications. A single-track program permits continued focus on a series of refereed submissions. A lively social program rounds out the event.

## Abstract

I make a naive implementation of the forward algorithm in Stan for the Normal Mixed GARCH. Using series for an index and stock prices from companies in different industries, I find that belief states are shared across assets and the strength of the relationship varies for each pair of assets. This hints that volatility states follow a hierarchical structure: for example, the risk states of a global portfolio may be decomposed in Country + Industry + Stock Individual.

## Foreword
The (very) **naive** implementation of the algorithm in Stan is only meant for illustration. A few good practices were neglected, convergence is not guaranteed, there is much room left for optimization and fitting *N* different independent models is probably not a reasonable choice for production sampler. The main takeaway of this presentation is the ideas behind the code but not the code itself.

## Prerequisites
  * R 3.3.3
  * RStudio Desktop 1.0.136
  * Rtools 3.3 (R 3.2.x to 3.3.x)
  * Stan 2.14
  * R Packages
    * RStan 2.14.2

## Authors

* **Luis Damiano** - [luisdamiano](https://github.com/luisdamiano)

## License

_A Quick Intro to Hidden Markov Models Applied to Stock Volatility_ is licensed under CC-BY-SA 4.0. See the [LICENSE](LICENSE.md) file for details.

## Acknowledgments

* To the R/Finance Conference committee for accepting my proposal and generously providing travel funding.
* Special thanks to all those who showed me how much fun stats can be, a real life changer.
