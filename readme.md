# Replication repository description for Nava and Dong (2022)'s JAAEA, *The impact of taxing sugar-sweetened beverages in México: A censored QUAI demand system approach*

Correspinding author: [Noé J Nava](noe.nava@usda.gov).

**Goal:** Replicate the results in `Table C.1`, `Table 2`, and `Table 3` in the paper. 

*Note:* Ou analysis is composed of two main procedues. A series of truncated log-likelihood functions are parametized and optimized to estimate the parameters described in `Table C.1`. Then, numerical approaches are employed, along witht the previously estimated parameters, to calculate demand elasticities shown in `Table 2` and `Table 3`. `Table 4` is superfluous to our main analysis, but such elasticities can be estimated using the same procedures employed for the previous tables.

<u>**`Table C.1`**:</u> There are two approaches to estimate the parameters in `Equation 9` and `Equation 12` as shown in `Table C.1` and described in `section 3.3` and `Appendix D`. The first is a tested one which reflects the results in our paper. The [APTECH: Gauss](https://www.aptech.com/) script is `code/TabC1_tested.gss` and requires the datasets `data/ssb_dataset_2018.dat` and `data/ssb_dataset_2018.dht`. The datasets are in [APTECH: Gauss](https://www.aptech.com/) format.

The second approach is *untested*. What it means is that we cannot guarantee the approach replicates the results in our paper since R has computational limitations. The script is `code/TabC1_untested.R` and uses the dataset `data/ssb_dataset_2018.csv`. The computational limitation is related to the estimation of the (n x p) Jacobian matrix required for our algorithm to find the solution. We were unable to check if our R script had any bug since our algorithm never converged.

<u>**`Table 2` and `Table 3`**:</u> Our numerical approach to the calculation of our demand elasticities is described in `Appendix E`. The R script is `code/Tab2_Tab3_tested.R` and uses output from the estimation of `code/TabC1_tested.gss`. Because I expect complications regarding the estimation of the parameters and the variance-covariance matrix depicted in `Table C.1`, the script is directed to use `TabC1_params.Rdata` for the parameters and `TabC1_vcov.xlsx` for the variance-covariance matrix.