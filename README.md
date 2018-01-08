# [Statistical anaylsis] New particle search in LHC big data
<div style="text-align: center;" markdown="1"><a href="http://dx.doi.org/10.1103/PhysRevD.93.112009">Phys. Rev. D 93 (2016) 112009</a></div>

*"Search for pair-produced vector-like quarks of charge -1/3 decaying to bH using boosted Higgs jet-tagging in pp collisions at sqrt(s) = 8 TeV"*

###### `big data` `background estimation` `statistics` `upper limits` `data-driven validation` `collabration` `particle physics`

## Introduction

The data anaylsis is for searching new pair-producted particles from proton-proton collision data collected by the CMS detector at CERN-LHC. The new particle is called vector-like quark (right- and left-chirality symmetry) with bottom flavor (**b'**), which is expected to be a possible mechanism to explaine the difference mass of the fundamental particles. Since the several data analysis results show the possible mass is out of the 800 GeV, this research focus on the boosted enegy scale.

In this research we focus on the decay of it with the bottom quark (**b**), which clusters a **jet** in data, and Higgs boson (**H**). Since the mass **b'** is predicted to be high, the heavy paritcle decay become possible. In contrast, the event (data) with heavy **b'** is extemly rare, e.g. one of millions events. The challange of this data analysis will be dominated by the background noice. Dealing with the **background estimation** and **validation** play the important roles.

On the other hand, when the mass of **b'** push to much higher scale, e.g. TeV, the **H** is expected to be boosted and dominantly decaying to a pair of bottom quarks ( the high decay probability is validated in both theory and experiment ). The pair of **b** from **H** will merge a boosted jets with a big cone shape, we called it **H jet**. Thus, the another important feature of this analysis is the feature extraction (indentification) of the boosted **H jet**.      

The overall new physics phenomenon is illustrated the following Figure:
<img src="http://hep1.phys.ntu.edu.tw/~alpha78718/cv/bprimetobH.png" height="300" >

## Techniques/results

#### 1. Data selection & feature selection
The data are recorded by several layer of detector and stored with electronic digits. The quality of data depends on the data taking processes and the reconstruction algorithm. The data must be filtered accodring to the condition of each processes, e.g. data taking trigger, performace vertexing algorithm of particle tracks etc....  

The features of the reconsturcted particles are also selected accodring to our anaylsis interests, e.g. at least 4 **jets** in data including at least 2 **b jet** and 1 high enegy **H jet**. Of course, the structure of dectors are also included, e.g. the particle is not possible to be detector from the gap of dectors. The follwing Figures show the predicted signal data and background noices, by ***Monte Carol (MC) method***, and data after selections. The backgound noices are expected to dominated in data by comparing the ratio beteew predicted signal data and background noice.

<img src="http://hep1.phys.ntu.edu.tw/~alpha78718/cv/bpbkg.png" height="200">

#### 2. Background estimation & validation
The background is expected to dominate in data from ***MC*** study. However, the background is actually hard to simulate, we have to estimate it from data. The used technique here for **predicting the background distrubtion (model)** of background noices in signal data is called ***Data-drven ABCD method***, which inverses the feature selections and obtains the model from the outlier data. Here I illustrated the method in this analysis, the signal region is at *B*, the backgound noices can be extroplated from *A* by applyig scale of *D/C*.

<img src="http://hep1.phys.ntu.edu.tw/~alpha78718/cv/bpabcd.png" height="200">  

The assumption of this method is expected the probability of backgorund noice is random and smoothly distrubted in overall data, and it has been validated by **MC method** and **Data-driven control sample** in this analysis. The contrl sample is also defined from the extremly outlier data, i.e. the data is not in *ABCD* regions. The selection of control sample requires the jets in data are not from **b**.

The total background in following Figure shows the predicted background model, which is matched with data, i.e. there is no any new significant new physics event (data) in current data-driven research.

<img src="http://hep1.phys.ntu.edu.tw/~alpha78718/cv/bpht1.png" height="250"> <img src="http://hep1.phys.ntu.edu.tw/~alpha78718/cv/bpht2.png" height="250">  

#### 3. Results
Although there is no any singal data found, it doesn't mean the physics mechanism is failed. We can still estimated the possible energy region which the new physics may exist there, i.e. the energy of LHC is not large enough to produce the new particle. Thus, we provide the **up-limit** of the probabilities of **b** piar production as function of mass of **b'** to shows the possible region to have **b'**.

The **up-limit** is estimated by [*Hypothesis Testing*](https://onlinecourses.science.psu.edu/statprogram/node/138) which give the conferdence interval with 95% by using data, predicted backgroud (from the ***Data-drven ABCD method***) and **b'**-MC samples. The following Figure shows the up-limit of the **b'** mass.

<img src="http://hep1.phys.ntu.edu.tw/~alpha78718/cv/bplimit.png" height="250">

- The blue dashed line shows the theoritial probabilities of **b** piar production which are applying to generate signal MC data.  

- The blue dot line within yellow and green bands present the predicted probabilities by data-drven predicted background and **b'**-MC sample. It compares the gaussian distrubtion of the yeilds of predicted background and prdiected background + **b'**-MC, i.e. comparing the predicted yeilds with and without singal exists. The probabilities of **b'** thus can be obtained by some caculation accoding to the coverage between two distrubtions. The yellow and green band are also estimated in same way but by variating the coverage. They are called *expected limit*.

- The black solid line is made by same way as expected limit, but the predicted background is replaced by data. This is called *observed limit*. Since there is no significant evidence excessing signal, the observe limit agrees with the expected limit within the uncertainty.  

- The widths of the qaussian distributions of the yields are obtained by uncertainty study. Except the statistical uncertainty, the systematic uncertainties are included due to the MC sample is used, and it is applied some corrections for matching with data. The detail study can be found in the **[References]**.

By comparing the probability of theory and limits, the observed (expectied) region of the **b'** mass can excludes the region below 846 (811) GeV at 95% confidence level.


## References
1. Github : <https://github.com/juifa-tsai/BprimeTobHAnalysis>
2. [Phys. Rev. D 93 (2016) 112009](http://dx.doi.org/10.1103/PhysRevD.93.112009).
3. [CERN-CMS-B2G-14-001](http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/B2G-14-001/index.html).
4. [CERN-CMS-B2G-13-006](http://cms-results.web.cern.ch/cms-results/public-results/publications/B2G-13-006/index.html)
5. [121th LHCC-CERN](https://indico.cern.ch/event/369822/#21-search-for-pair-produced-ve)
