---
title: "Towards combined epidemiological and economic dynamical systems"
#date: "24 February, 2025"
output:
  bookdown::github_document2:
    # pandoc_args: --webtex
    toc: true
    toc_depth: 5
    toc_float: true
    number_sections: true
  bookdown::pdf_document2: 
    toc: false
    keep_tex: yes
    citation_package: natbib
    extra_dependencies: ["float"]
    # extra_dependencies: ["flafter"]
    pandoc_args:
      --filter=pandoc-xnos
    number_sections: true
    fig_caption: yes
    includes:
      # in_header: "preamble.tex"
  bookdown::word_document2: 
    toc_depth: 5
    toc_float: true
    number_sections: true
    editor_options: 
      chunk_output_type: inline
bibliography: 
  - "DAEDALUS.bib"
header-includes:
always_allow_html: true
---

This document describes some combined epi-econ models, starting with
DAEDALUS, and then moving onto SFC-epi models.

Code examples are taken and adapted from
https://github.com/marcoverpas/Italy-SFC-Model and
https://github.com/marcoverpas/Six_lectures_on_sfc_models.

# [1 DAEDALUS]{data-rmarkdown-temporarily-recorded-id="daedalus"}

Originally designed of for the UK [@Haw2020] and since applied to
Indonesia [@Johnson2023]. The econ model is static. The industrial
sectors of the economy are used to stratify the population in the epi
model, in which there are contacts related to sector of work and
consumption.

In response to an outbreak, closures are mandated, resulting in
operations in each sector reducing to some percentage. In the model, GVA
reduces to the same percentage, and so does attendance of places of work
and consumption, which scales contacts made and subsequently disease
transmission.

The population is static: it does not grow, and no one changes sector of
work. The non-working group remains non-working. The workforce in place
includes a fraction who work from home and therefore contribute to GVA
but not infections from workplace interactions.

See https://github.com/robj411/p2_drivers for a recent presentation.

::: figure
`<img src="README_files/figure-gfm/dae-1.png" alt="Daedalus model. Mandated closures reduce both infections and workforce in place." width="70%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:dae}Figure 1.1: Daedalus model. Mandated closures reduce both
infections and workforce in place.
```{=html}
</p>
```
:::

# [2 SFC-epi]{data-rmarkdown-temporarily-recorded-id="sfc-epi"}

We will use an SFC model in place of GVA by sector. The SFC model is
dynamic in time, meaning impacts can accumulate, and is demand driven,
meaning we can model autonomous behaviour changes that will have impacts
on both the economy and the epidemic.

At present we compute the "propensity to consume/work" which scales the
contacts in the epi model and feeds into the econ model. It might be
better to precompute the counterfactual (without epidemic) consumption
and labour, and:

1.  compute scaled parameters using epi variable(s) and response
    function;
2.  solve econ model to get consumption and labour;
3.  scale contacts with the same ratio as consumption and labour have to
    the no-epidemic counterfactual at the same time point

NB: the model assumes an impact on the general population (i.e. both
susceptibles and infectious change their behaviour). Lack of consumption
and work from people who are already sick *because* they are sick is not
modelled. This could be included explicitly by e.g. setting the
proportion of non-diseased people as the upper bound to propensities.

::: figure
`<img src="figures/response.png" alt="Dose--response function: extent of behaviour change as a function of an epidemiological variable such as number of hospital cases." width="50%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:unnamed-chunk-1}Figure 2.1: Dose--response function: extent of
behaviour change as a function of an epidemiological variable such as
number of hospital cases.
```{=html}
</p>
```
:::

## [2.1 Consumption]{data-rmarkdown-temporarily-recorded-id="consumption"}

::: figure
`<img src="README_files/figure-gfm/consumption-1.png" alt="Epi variables reduce propensity to consume, which reduces consumption, and therefore GVA and new infections." width="70%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:consumption}Figure 2.2: Epi variables reduce propensity to
consume, which reduces consumption, and therefore GVA and new
infections.
```{=html}
</p>
```
:::

The first step is to model consumption as a function of the epidemic,
e.g. the number of cases, hospitalisations or deaths. The epidemic
variable modifies the "propensity to consume" parameters (often labelled
$\alpha$). We scale propensity to consume and exposure-related
activities by the same amount.

::: figure
`<img src="figures/SIM_1e+05-0.5.png" alt="Results for SIM model with consumption reduction alongside counterfactual epi and econ curves without integration." width="80%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:unnamed-chunk-2}Figure 2.3: Results for SIM model with
consumption reduction alongside counterfactual epi and econ curves
without integration.
```{=html}
</p>
```
:::

## [2.2 Labour supply]{data-rmarkdown-temporarily-recorded-id="labour-supply"}

::: figure
`<img src="README_files/figure-gfm/labour-1.png" alt="Epi variables reduce propensity to work and to consume, which reduces consumption and labour supply, and therefore GVA and new infections." width="70%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:labour}Figure 2.4: Epi variables reduce propensity to work and
to consume, which reduces consumption and labour supply, and therefore
GVA and new infections.
```{=html}
</p>
```
:::

Next, we model reduction in propensity to work in the same way.

However, this should ultimately depend also on the need for income,
meaning that (a) lower income people are more likely to work, (b) people
are more likely to work as time goes on and they've spent some months
with less income, and (c) we can later mitigate these effects to some
extent with government transfers.

This version of the model (with spontaneous consumption and labour
change) can be used to model an unmitigated epidemic: that is, an
epidemic where there are no government interventions, but the population
can choose infection-avoiding behaviours, which will dampen both the
epidemic and the economy. This is an important gap in current epi
modelling. We usually assume either no behaviour change, so that
unmitigated epidemics have huge death tolls (e.g. UK COVID projection),
or that uncosted (free) population behaviour change will curb the
excesses of the epidemic (CEPI work).

## [2.3 Mandated closures]{data-rmarkdown-temporarily-recorded-id="mandated-closures"}

::: figure
`<img src="README_files/figure-gfm/mandate-1.png" alt="Epi variables reduce propensity to work and to consume, which reduces consumption and labour supply, and therefore GVA and new infections. Mandated closures reduce consumption and propensity to work." width="70%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:mandate}Figure 2.5: Epi variables reduce propensity to work and
to consume, which reduces consumption and labour supply, and therefore
GVA and new infections. Mandated closures reduce consumption and
propensity to work.
```{=html}
</p>
```
:::

Next we want to introduce the ability for the government to mandate
closures. This might be modelled as e.g. a maximum labour demand,
expressed as a percentage of counterfactual labour demand.

Propose a worked example, such as "how would we model a mandate that the
hospitality sector closes?"

## [2.4 Government transfers]{data-rmarkdown-temporarily-recorded-id="government-transfers"}

::: figure
`<img src="README_files/figure-gfm/transfers-1.png" alt="Epi variables reduce propensity to work and to consume, which reduces consumption and labour supply, and therefore GVA and new infections. Mandated closures reduce consumption and propensity to work. Government transfers increase consumption and reduce propensity to work." width="70%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:transfers}Figure 2.6: Epi variables reduce propensity to work
and to consume, which reduces consumption and labour supply, and
therefore GVA and new infections. Mandated closures reduce consumption
and propensity to work. Government transfers increase consumption and
reduce propensity to work.
```{=html}
</p>
```
:::

We should introduce government transfers to demonstrate that mandated
closures are only sustainable for as long as the population is supported
to forego income in order to stop the spread of infection.

## [2.5 Other things to consider]{data-rmarkdown-temporarily-recorded-id="other-things-to-consider"}

-   international trade, esp. tourism
-   structural changes over time?
-   move to online consumption

# [3 Epidemic model]{data-rmarkdown-temporarily-recorded-id="epidemic-model"}

The epidemic model is similar to DAEDALUS:

-   four age groups (pre-school age, school-age children, working age,
    retirement age)
-   working-age people stratified by something: definitely working/not
    working, potentially something sector, occupation or SES related
-   seven disease states (S, E, I (symptomatic/asymptomatic), R, H, D)
-   no vaccination or waning, for simplicity
-   leave school closures out for now, for simplicity
-   outputs from the econ model parametrise the function $k_{j}^{1}$

## [3.1 Ordinary differential equations]{data-rmarkdown-temporarily-recorded-id="ordinary-differential-equations"}

$$\begin{align}
\frac{dS_{j}}{dt} & = - k_{j}^{1}(t)S_{j}  \\
\frac{dE_{j}}{dt} & = k_{j}^{1}(t)S_{j} - (k^2+k^4)E_{j} \\
\frac{dI_{j}^a}{dt} & = k^2E_{j} - k^3I_{j}^a \\
\frac{dI_{j}^s}{dt} & = k^4E_{j} - (k_{j}^{5}+k_{j}^{6})I_{j}^s \\
\frac{dR_{j}}{dt} & = k^3I_{j}^a + k_{j}^{5}I_{j}^s + k_{j}^{7}(t) H_{j}\\
\frac{dH_{j}}{dt} & = k_{j}^{6}I_{j}^s - (k_{j}^{7}(t) + k_{j}^{8}(t)) H_{j} \\
\frac{dD_{j}}{dt} & =  k_{j}^{8}(t) H_{j}
\end{align}$$

## [3.2 Disease state transitions]{data-rmarkdown-temporarily-recorded-id="disease-state-transitions"}

::: figure
`<img src="README_files/figure-gfm/statetransitions-1.png" alt="Disease state transitions. $S$: susceptible. $E$: exposed. $I^{a}$: asymptomatic infectious. $I^{s}$: symptomatic infectious. $H$: in need of hospitalisation. $R$: recovered. $D$: died. $j$: stratum." width="50%" />`{=html}
```{=html}
<p class="caption">
```
[]{#fig:statetransitions}Figure 3.1: Disease state transitions. $S$:
susceptible. $E$: exposed. $I^{a}$: asymptomatic infectious. $I^{s}$:
symptomatic infectious. $H$: in need of hospitalisation. $R$: recovered.
$D$: died. $j$: stratum.
```{=html}
</p>
```
:::

# [4 References]{data-rmarkdown-temporarily-recorded-id="references"}
