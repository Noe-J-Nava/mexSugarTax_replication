# Replication repository description for Nava and Dong (2022)'s JAAEA, *The impact of taxing sugar-sweetened beverages in México: A censored QUAI demand system approach*

Correspinding author: [Noé J Nava](noe.nava@usda.gov).

**Goal:** Replicate the results in `Table C.1`, `Table 2`, and `Table 3` in the paper. 

*Note:* Our analysis is composed of two main procedures. A series of truncated log-likelihood functions are parametized and optimized to estimate the parameters described in `Table C.1`. Then, numerical approaches are employed, along with the previously estimated parameters, to calculate demand elasticities shown in `Table 2` and `Table 3`. `Table 4` is superfluous to our main analysis, but such elasticities can be estimated using the same procedures employed for the previous tables and previous data.

### `Table C.1`: 

There are two approaches to estimate the parameters in `Equation 9` and `Equation 12` as shown in `Table C.1` and described in `section 3.3` and `Appendix D`. The first is a tested one which reflects the results in our paper. The [APTECH: Gauss](https://www.aptech.com/) script is `code/TabC1_tested.gss` and requires the datasets `data/ssb_dataset_2018.dat` and `data/ssb_dataset_2018.dht`. The datasets are in [APTECH: Gauss](https://www.aptech.com/) format.

The second approach is *untested*. What it means is that we cannot guarantee the approach replicates the results in our paper since R has computational limitations. The script is `code/TabC1_untested.R` and uses the dataset `data/ssb_dataset_2018.csv`. The computational limitation is related to the estimation of the (n x p) Jacobian matrix required for our algorithm to find the solution. We were unable to check if our R script had any bug since our algorithm never converged. Despite that, we prepared this R script for those interested in the computational task behind our paper but do not know how to read Gauss scripts.

I am testing Julia language at the moment, but I was told the language can address the computational task.

### `Table 2` and `Table 3`:

Our numerical approach to the calculation of our demand elasticities is described in `Appendix E`. The R script is `code/Tab2_Tab3_tested.R` and uses output from the estimation of `code/TabC1_tested.gss`. Because I expect complications regarding the estimation of the parameters and the variance-covariance matrix depicted in `Table C.1`, the script is directed to use `TabC1_params.Rdata` for the parameters and `TabC1_vcov.xlsx` for the variance-covariance matrix.

**How to cite this article:**

Nava, Noé J., and Diansheng Dong. 2022. "The impact of taxing sugary-sweetened beverages in México: A censored QUAI demand system approach." *Journal of Agricultural and Applied Economics Association*, 1(1):1-23: https://doi.org/10.1002/jaa2.6

