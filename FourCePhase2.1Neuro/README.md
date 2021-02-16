Analysis of COVID-19 patients with Neurological Conditions
================

## Install and run

First, on your command line, open the Docker container, replacing
`your_path_here` with the path to the directory containing folder
`Input` (where your site-specific datasets are stored):

``` bash
docker run \
  --name 4ce \
  --volume your_path_here:/4ceData \
  --rm -it \
  dbmi/4ce-analysis:version-2.0.0 R
```

This should open up the interactive R environment on your terminal. Now,
while in R, install the development version of **Phase2.1NeuroRPackage**
from GitHub with remotes:

``` r
# install.packages('remotes') # uncomment to install devtools
remotes::install_github('covidclinical/Phase2.1NeuroRPackage',
                        subdir = 'FourCePhase2.1Neuro',
                        upgrade = FALSE)
```

The main function `runAnalysis()` has 4 required arguments. See your
site obfuscation parameters (`mask_thres` and `blur_abs`)
[here](https://docs.google.com/spreadsheets/d/1Xl9juDBXt86P3xQtsoTaBl2zPl1BIiAG9DI3Rotyqp8/edit#gid=212461777).
If the site’s ICD codes are primarily of version 9, set
`icd_version = 9`. If your site does not have race information, set
`include_race = FALSE`.

``` r
library(Phase2.1NeuroRPackage)

runAnalysis(mask_thres = 10, blur_abs = 0, icd_version = 10, include_race = TRUE)
```

Finally, please submit the results to the central repository. Please
note that you would need to use a *token* to access private repos, see
[here](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token).

``` r
submitAnalysis()
```

Briefly, to generate a new token, go to your GitHub settings -&gt;
Developer settings -&gt; Personal access tokens -&gt; Generate.

![](images/token.png)

## Other notes

If everything runs smoothly, your result would be under `Phase2.1NeuroR`
in your home directory, *i.e.*, `~/Phase2.1NeuroR`. For example, for me,
it’s `/data2/home/ttle/Phase2.1NeuroR/penn_results.rda`. If somehow
`submitAnalysis()` didn’t allow you to upload the results to
[Phase2.1NeuroRSummariesPublic](https://github.com/covidclinical/Phase2.1NeuroRSummariesPublic)
, you can share the results file with us (@meghutch and @trang1618) via
the \#neuro Slack channel.

To get back to the command line at any point, type `quit()` in the R
environment. Note that this will exit the container, and you would have
to open another one and reinstall the package if you want to rerun.

If you run into any problem adapting this code to your data, let us know
via Slack or [submit an
issue](https://github.com/covidclinical/Phase2.1NeuroRPackage/issues/new).

Thank you very much for your support!
