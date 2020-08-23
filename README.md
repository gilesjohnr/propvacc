[![Build Status](https://travis-ci.com/gilesjohnr/propvacc.svg?branch=master)](https://travis-ci.com/gilesjohnr/propvacc)

---
title: "General theory to estimate doses of vaccine received and proportion of population immune due to vaccination"
author: "John Giles"
date: '15 May 2020'
output:
  html_document: default
  pdf_document: default
header-includes:
- \usepackage{amsmath}
- \usepackage{amssymb}
- \usepackage{amsfonts}
editor_options: 
  chunk_output_type: inline
---

This repository contains the propvacc R package which provides functions to model and simulate vaccine dose coverage.

# Introduction
Vacine preventable diseases are important since we have highly effacacious vaccines for things like Polio, measles, ruebella etc, but inefficient distribution of these leads to gaps in susceptibility that result in outbreaks. One major inefficiency is re-vaccinating children so that some have more than the required number of doses, but the unvaccinated population remains untouched. Optimal distribution requires understanding where unvaccinated or under-vaccinated individuals remain (in particular those that have received zero doses). Targeting those with zero-doses is vital countries approaching elimination.

This is particularly challenging because our current methods often over-estimate the proportion of the population that is immune due to vaccination. Sero-surveys are an accurate way to measure overall immunity, however, these include recovered individuals and cannot be used to estimate proportion vaccinated in a straightforward manner. Modeling studies have used DHS and WHO data to estimate the proportion vaccinated and often assume independence between each dose received from routine vaccination and/or SIA campaigns. However, in practice vaccine doses administered during subsequent routine activities are often dependent on whether previous doses were received. There are studies that try to quantify the extent of dependence on prior vaccination, however, they are few and limited to countries that have performed post campiagn surveys and added questions regarding prior vaccination. Therefore a major challenge is modeling the dependence among routine activities and campaigns when only independent measures of coverage are available.

Here, we develop a parameter-free model of population immunity that uses independent measurements of vaccination coverage and accounts for dependence among doses. This method also estimates the proportion of the population that has received each possible number of doses given vaccination activities by defining conditional probabilities based on a dropout rate between doses and using the law of total probability to calculate relative proportions of each dose. For example, in the case of two-dose routine vaccination, we assume that the probability of receiving dose two regardless of age is higher if an individual already received dose one. We then estimate the proportion of the population that has received 0, 1, and 2 doses of the vaccine.

# Methods
To use general terms, we define the proportion successfully vaccinated in each of the $k$ vaccination activities using the vector $\{ V_1, V_2, \cdots, V_k \}$ and $S$ is the proportion vaccinated in the SIA campaign. 

## Estimating proportion vaccinated for a two-dose vaccine with an SIA campaign
Assumption of dependence:

  1. Successful vaccination in routine activity is dependent upon the number of doses received previously, where the coverage of the second vaccination activity is first distributed to those that have one prior dose and then the surplus (if any) is distributed to those that have zero prior doses.
  2. Vaccination in an SIA campaign is higher for those that have received at least one prior dose ([Subaiya et al 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6028100/pdf/pone.0199786.pdf)). Where, the coverage of the SIA campaign is first distributed to those with at least one prior dose and then the surplus (if any) is distributed to those that have zero prior doses.

### Proportion that have received at least one dose prior to SIA campaign
The expected proportion of the population that has received one or two doses of vaccine prior to the SIA campaign is:
$$
\begin{aligned}
\Pr(\text{prior dose})
&=
\frac{ 
\Pr( V_1 , V_2) + \Pr( V_1 , \neg V_2 ) + \Pr( \neg V_1 , V_2)
}{  
\Pr( V_1 , V_2) + \Pr( V_1 , \neg V_2 ) + \Pr( \neg V_1 , V_2) + \Pr( \neg V_1 , \neg V_2)
}
\\[10pt]
\end{aligned}
$$

The conditional probability terms are defined as:
$$
\begin{aligned}
\Pr( V_1 , V_2)  &= \text{probability of receiving dose one and dose two}\\[5pt]
\Pr( V_1 , \neg V_2 ) &= \text{probability of receiving dose one and missing dose two}\\[5pt]
\Pr( \neg V_1, V_2) &= \text{probability of missing dose one and receiving dose two}\\[5pt]
\Pr( \neg V_1 , \neg V_2) &= \text{probability of missing  both dose one and dose two}\\[5pt]
\end{aligned}
$$

If we assume that individuals that received dose one are the most likely to receive dose two, we can define the dropout rate as in ([Masresha et al 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6372060/)) to get the conditional probability of receiving dose two after dose one: $d_{1,2}$. Since the proportion of dose two can technically be higher than dose one (albeit unlikely), the dropout out term becomes zero when $V_2 \ge V_1$:

$$
d_{1,2} = \begin{cases}
          0, & \text{if} \ V_2 \ge V_1. \\[5pt]
          \frac{(V_1 - V_2)}{V_1}, & \text{if} \ V_2 < V_1.
          \end{cases}
$$
Under this particular assumption of dependence, the conditional probability terms used to calculate the probability of receiving any doses prior to an SIA campaign can be calculated as:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( V_1 , V_2) \big] = \Pr( V_2 \mid V_1) \Pr( V_1 ) &= V_1(1-d_{1,2}) 
\\[5pt]
\mathbb{E}\big[\Pr( V_1 , \neg V_2) \big] = \Pr( \neg V_2 \mid V_1) \Pr( V_1 ) &= V_1d_{1,2} 
\\[5pt]
\mathbb{E}\big[\Pr( \neg V_1 , V_2) \big]= \Pr( V_2 \mid \neg V_1) \Pr( \neg V_1 ) &= 
          \begin{cases}
          0, & \text{if} \ V_2 \leq V_1. \\
          (1-V_1)(V_2-V_1), & \text{if} \ V_2 > V_1.
          \end{cases}
\\[5pt]
\mathbb{E}\big[\Pr( \neg V_1 , \neg V_2) \big]= \Pr( \neg V_2 \mid \neg V_1) \Pr( \neg V_1 ) &= 
          \begin{cases}
          1-V_1, & \text{if} \ V_2 \leq V_1. \\
          (1-V_1)(1-(V_2-V_1) ), & \text{if} \ V_2 > V_1.
          \end{cases}
\end{aligned}
$$

### Incorporating SIA campaign into two-dose model
Building on the above assumption of dependence between $V_1$ and $V_2$, the conditional probability of each dose based on two-dose routine vaccination acitivity, we add a term for the SIA campaign $S$, which is dependent on whether at least one prior dose was received. The denominator is comprised of all potential outcomes of two routine vaccination activities and one SIA campaign. The total sample space of the two-dose model and SIA campaign is comprised of 8 unique combinations, which I've written as $\Omega_{V_1, V_2, S}$:
$$
\Omega_{V_1, V_2, S} = 
\begin{bmatrix}
\Pr( V_1 , V_2, S) + \Pr( V_1 , \neg V_2, S ) + \Pr( \neg V_1 , V_2, S) + \Pr( \neg V_1 , \neg V_2, S) +\ \\
\Pr( V_1 , V_2, \neg S) + \Pr( V_1 , \neg V_2, \neg S ) + \Pr( \neg V_1 , V_2, \neg S) + \Pr( \neg V_1 , \neg V_2, \neg S)
\end{bmatrix}
$$
And the probability terms for each dose are:
$$
\begin{aligned}
\Pr(\text{three doses})
&=
\frac{ 
\Pr( V_1 , V_2, S)
}{  
\Omega_{V_1, V_2, S}
}
\\[10pt]
\Pr(\text{two dose})
&=
\frac{ 
\Pr( \neg V_1 , V_2, S) + \Pr( V_1 , \neg V_2, S ) + \Pr( V_1 , V_2, \neg S)
}{  
\Omega_{V_1, V_2, S}
}
\\[10pt]
\Pr(\text{one dose})
&=
\frac{ 
\Pr( \neg V_1 , \neg V_2, S) + \Pr( V_1 , \neg V_2, \neg S ) + \Pr( \neg V_1 , V_2, \neg S)
}{  
\Omega_{V_1, V_2, S}
}
\\[10pt]
\Pr(\text{zero doses})
&=
\frac{ 
\Pr( \neg V_1 , \neg V_2, \neg S)
}{  
\Omega_{V_1, V_2, S}
}
\end{aligned}
$$

Using our definition of $\Pr(\text{prior dose})$ above we can calculate the probability of receiving a dose through an SIA campaign conditionaed on prior vaccination by defining the dropout rate from prior vaccination to SIA campaign:
$$
d_S = \begin{cases}
          0, & \text{if} \ S \ge \Pr(\text{prior dose}). \\[5pt]
          \frac{(\Pr(\text{prior dose}) - S)}{\Pr(\text{prior dose})}, & \text{if} \ S < \Pr(\text{prior dose}).
          \end{cases}
$$

Under this particular assumption of dependence between routine activity and SIA campaign coverage, the conditional probability terms defined above can then be calculated as:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( V_1, V_2, S) \big] = \Pr( V_2 \mid V_1) \Pr( V_1 ) \Pr( S \mid \text{prior dose} ) &= V_1(1-d_{1,2})(1-d_S)
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , V_2, \neg S) \big] = \Pr( V_2 \mid V_1) \Pr( V_1 ) \Pr( \neg S \mid \text{prior dose}) &= V_1(1-d_{1,2})d_S
\\[5pt]
\mathbb{E}\big[\Pr( V_1 , \neg V_2, S) \big] = \Pr( \neg V_2 \mid V_1) \Pr( V_1 ) \Pr( S \mid \text{prior dose}) &= V_1d_{1,2}(1-d_S)
\\[5pt]
\mathbb{E}\big[\Pr( V_1 , \neg V_2, \neg S) \big] = \Pr( \neg V_2 \mid V_1) \Pr( V_1 ) \Pr(\neg S \mid \text{prior dose}) &= V_1d_{1,2}d_S 
\\[5pt]
\mathbb{E}\big[\Pr( \neg V_1 , V_2, S) \big]= \Pr( V_2 \mid \neg V_1) \Pr( \neg V_1 ) \Pr(S \mid \text{prior dose}) &= 
          \begin{cases}
          0, & \text{if} \ V_2 \leq V_1. \\[3pt]
          (1-V_1)(V_2-V_1)(1-d_S), & \text{if} \ V_2 > V_1.
          \end{cases}
\\[5pt]
\mathbb{E}\big[\Pr( \neg V_1 , V_2, \neg S) \big]= \Pr( V_2 \mid \neg V_1) \Pr( \neg V_1 ) \Pr(\neg S \mid \text{prior dose}) &= 
          \begin{cases}
          0, & \text{if} \ V_2 \leq V_1. \\[3pt]
          (1-V_1)(V_2-V_1)d_S, & \text{if} \ V_2 > V_1.
          \end{cases}
\\[5pt]
\mathbb{E}\big[\Pr( \neg V_1 , \neg V_2, S) \big]= \Pr( \neg V_2 \mid \neg V_1) \Pr( \neg V_1 ) \Pr(S \mid \text{prior dose}) &= 
          \begin{cases}
          (1-V_1)(1-d_S), & \text{if} \ V_2 \leq V_1. \\[3pt]
          (1-V_1)(1-(V_2-V_1) )(1-d_S), & \text{if} \ V_2 > V_1.
          \end{cases}\\[5pt]
\mathbb{E}\big[\Pr( \neg V_1 , \neg V_2, \neg S) \big] = \Pr( \neg V_2 \mid \neg V_1) \Pr( \neg V_1 ) \Pr(\neg S \mid \text{prior dose}) &= 
          \begin{cases}
          (1-V_1)d_S, & \text{if} \ V_2 \leq V_1. \\[3pt]
          (1-V_1)(1-(V_2-V_1) )d_S, & \text{if} \ V_2 > V_1.
          \end{cases}
\end{aligned}
$$

### Total proportion immune due to vaccination
Given the definitions of recieving 0, 1, or 2 doses shown above, the total proportion vaccinated $p_{\text{vacc}}$ including effectiveness is:
$$
p_{\text{vacc}} = \sum_j \text{effectiveness}_j \times \Pr(\text{dose}_j), \ \ \text{where} \  j = \{ 1, 2, 3 \}
$$



## Estimating proportion vaccinated for a three-dose vaccine with an SIA campaign
Assumption of dependence:

  1. Successful vaccination in routine activity is dependent upon the number of doses received previously, where the coverage of the third vaccination activity is first distributed to those that have two prior doses, and then the surplus (if any) is distributed to those that have one prior dose, and again, if there is still a surplus, it is distributed to those with zero prior doses.
  2. Vaccination in the SIA campaign is the same as the two-dose method, where it is higher for those that have received one or more prior doses ([Subaiya et al 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6028100/pdf/pone.0199786.pdf)). Where, the coverage of the SIA campaign is first distributed to those with any number of prior doses and then the surplus (if any) is distributed to those that have zero prior doses.

#### Proportion that have received any doses prior to SIA campaign
Here we estimate the proportion of the population that has received 0, 1, 2, or 3 doses of a vaccine. We define the proportion successfully vaccinated with each round of vaccination as $\{ V_1, V_2, V_3\}$. Where, the probability of vaccination with dose $V_k$ is estimated using independent observations of the proportion of the population that has received dose $k$.


### Proportion that have received at least one dose prior to SIA campaign
The total sample space of a routine activities of a three dose vaccine is comprised of 8 unique combinations which I've written as $\Omega_{V_1, V_2, V_3}$:
$$
\Omega_{V_1, V_2, V_3} = \begin{bmatrix}
\Pr( V_1 , V_2, V_3) + \Pr( V_1 , V_2, \neg V_3) + \Pr( V_1 , \neg V_2, V_3) + \Pr( \neg V_1 , V_2, V_3) +\ \\
\Pr( V_1 , \neg V_2, \neg V_3) + \Pr( \neg V_1 , V_2, \neg V_3) + \Pr( V_1 , \neg V_2, \neg V_3) + \Pr( \neg V_1 , \neg V_2, \neg V_3) \\
\end{bmatrix}
$$

The expected proportion of the population that has received at least one dose through routine activity prior to an SIA campaign is then:
$$
\begin{aligned}
\Pr(\text{prior dose})
&=
\frac{ 
\Omega_{V_1, V_2, V_3} - \Pr( \neg V_1 , \neg V_2, \neg V_3) 
}{  
\Omega_{V_1, V_2, V_3}
}
\end{aligned}
$$

The conditional probability terms are defined as:

$$
\begin{aligned}
\Pr( V_1 , V_2, V_3) &= \text{probability of receiving all three doses}\\[5pt]
\Pr( V_1 , V_2, \neg V_3) &= \text{probability of receiving doses one and two and missing dose three}\\[5pt] 
\Pr( V_1 , \neg V_2, V_3) &= \text{probability of receiving doses one and three and missing dose two}\\[5pt]
\Pr( \neg V_1 , V_2, V_3) &= \text{probability of receiving doses two and three and missing dose one}\\[5pt]
\Pr( V_1 , \neg V_2, \neg V_3) &= \text{probability of missing doses two and three and receiving dose one}\\[5pt]
\Pr( \neg V_1 , V_2, \neg V_3) &= \text{probability of missing doses one and three and receiving dose two}\\[5pt]
\Pr( \neg V_1 , \neg V_2, V_3) &= \text{probability of missing doses one and two and receiving dose three}\\[5pt] 
\Pr( \neg V_1 , \neg V_2, \neg V_3) &= \text{probability of missing all three doses}\\[5pt]
\end{aligned}
$$
Defining the dropout rate from dose 2 to dose 3:

If we assume that individuals that received doses 1 and 2 are the most likely to receive dose 3, we can define the dropout rate $d_{2,3}$ to get the conditional probability of receiving dose 3 after receiving both dose 1 and 2. Since the proportion of dose 3 can technically be higher than those that have receivied doses 1 and 2, the dropout out term becomes zero when $V_3 \ge \Pr(V_2, V_1)$:

$$
d_{2,3} = \begin{cases}
          0, & \text{if}  \ V_3 \geq \Pr(V_2, V_1). \\[3pt]
          \frac{( \Pr(V_2 \mid V_1) - V_3 )}{\Pr(V_2 \mid V_1)}, & \text{if} \ V_3 < \Pr(V_2, V_1).
          \end{cases}
$$

Under this particular assumption of dependence, the expected values of the probability terms are:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( V_1 , V_2, V_3) \big] = 
\Pr(V_3 \mid V_1,V_2) \Pr( V_2 \mid V_1) \Pr( V_1 ) &= V_1(1-d_{1,2})(1-d_{2,3}) 
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , V_2, \neg V_3) \big] = 
\Pr(\neg V_3 \mid V_1,V_2) \Pr( V_2 \mid V_1) \Pr( V_1 ) &= V_1(1-d_{1,2})d_{2,3} 
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , \neg V_2, V_3) \big] = 
\Pr(V_3 \mid V_1, \neg V_2) \Pr(\neg V_2 \mid V_1) \Pr(V_1) &= 
     \begin{cases}
     0,                                                & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\[3pt]
     V_1 d_{1,2} \big(V_3 - (V_1(1-d_{1,2})) \big),    & \text{if} \ V_3 > \Pr(V_2, V_1).
     \end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , V_2, V_3) \big] = 
\Pr(V_3 \mid \neg V_1, V_2) \Pr( V_2 \mid \neg V_1) \Pr(\neg V_1) &= 
     \begin{cases}
     0,                                                               & \text{if} \ V_2 \leq V_1. \\[3pt]
     0,                                                               & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\[3pt]
     \big((1-V_1)(V_2-V_1)\big)  \big(V_3 - (V_1(1-d_{1,2})) \big),   & \text{otherwise}.
     \end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , \neg V_2, \neg V_3) \big] =
\Pr(\neg V_3 \mid V_1, \neg V_2) \Pr( \neg V_2 \mid V_1) \Pr(V_1) &= 
     \begin{cases}
     V_1 d_{1,2},                                                & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\[3pt]
     V_1 d_{1,2} \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big),        & \text{if} \ V_3 >  \Pr(V_2, V_1). 
     \end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , V_2, \neg V_3) \big] =
\Pr(\neg V_3 \mid \neg V_1, V_2) \Pr( V_2 \mid \neg V_1) \Pr(\neg V_1) &= 
     \begin{cases}
     0,                                                               & \text{if} \ V_2 \leq V_1. \\[3pt]
     0,                                                               & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\[3pt]
     (1-V_1)(V_2-V_1) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big),        & \text{otherwise}.
     \end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , \neg V_2, V_3) \big] =
\Pr(V_3 \mid \neg V_1, \neg V_2) \Pr( \neg V_2 \mid \neg V_1) \Pr(\neg V_1) &= 
     \begin{cases}
      0,                                                                                                  & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\[3pt]
      (1-V_1) \big( V_3 - V_1d_{1,2} - (1-V_1)(V_2-V_1) - V_1(1-d_{1,2}) \big),                           & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\[3pt]
      (1-V_1) \big(1-(V_2-V_1)\big) \big( V_3 - V_1d_{1,2} - (1-V_1)(V_2-V_1) - V_1(1-d_{1,2}) \big),     & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 >  V_1.
     \end{cases}
\\[5pt] 
\mathbb{E}\big[ \Pr( \neg V_1 , \neg V_2, \neg V_3) \big] =
\Pr(\neg V_3 \mid \neg V_1, \neg V_2) \Pr( \neg V_2 \mid \neg V_1) \Pr(\neg V_1) &= 
     \begin{cases}
      1-V_1,                                                                    & \text{if} \ V_3 \leq \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\[3pt]
      (1-V_1) \big(1-(V_2-V_1)\big),                                            & \text{if} \ V_3 \leq \Pr(V_2, V_1) \ \& \ V_2 >  V_1. \\[3pt]
      (1-V_1) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big),                          & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\[3pt]
      (1-V_1) \big(1-(V_2-V_1)\big) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big),    & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 >  V_1.
     \end{cases}
\end{aligned}
$$


### Incorporating SIA campaign into three-dose model
The denominator is comprised of all potential outcomes of the three routine activities and one SIA campaign. The total sample space of the three-dose model and SIA campaign is comprised of 16 unique combinations, which I've written as $\Omega_{V_1, V_2, V_3, S}$:
$$
\Omega_{V_1, V_2, V_3, S} = \begin{bmatrix}
\Pr( V_1 , V_2, V_3, S) + \Pr( V_1 , V_2, \neg V_3, S) + \Pr( V_1 , \neg V_2, V_3, S) + \Pr( \neg V_1 , V_2, V_3, S) +\ \\
\Pr( V_1 , \neg V_2, \neg V_3, S) + \Pr( \neg V_1 , V_2, \neg V_3, S) + \Pr( V_1 , \neg V_2, \neg V_3, S) + \Pr( \neg V_1 , \neg V_2, \neg V_3, S) +\ \\
\Pr( V_1 , V_2. V_3, \neg S) + \Pr( V_1 , V_2, \neg V_3, \neg S) + \Pr( V_1 , \neg V_2, V_3, \neg S) + \Pr( \neg V_1 , V_2, V_3, \neg S)+\  \\
\Pr( V_1 , \neg V_2, \neg V_3, \neg S) + \Pr( \neg V_1 , V_2, \neg V_3, \neg S) + \Pr( V_1 , \neg V_2, \neg V_3, \neg S) + \Pr( \neg V_1 , \neg V_2, \neg V_3, \neg S) \\
\end{bmatrix}
$$

Using the above sample space $\Omega_{V_1, V_2, V_3, S}$, we can write the conditional probabilies for receiving 0, 1, 2, 3, or 4 doses by building on the assumptions of dependence among routine activies $V_1$, $V_2$, and $V_3$, and one SIA campaign $S$:
$$
\begin{aligned}
\Pr(\text{four doses}) &=
\frac{ 
\Pr( V_1, V_2, V_3, S)
}{  
\Omega_{V_1, V_2, V_3, S}
}
\\[10pt]
\Pr(\text{three doses})
&=
\frac{ 
\Pr(V_1, V_2, V_3, \neg S) + \Pr(V_1 , V_2, \neg V_3, S) + \Pr(V_1 , \neg V_2, V_3, S) + \Pr(\neg V_1 , V_2, V_3, S) 
}{  
\Omega_{V_1, V_2, V_3, S}
}
\\[10pt]
\Pr(\text{two doses})
&=
\frac{ 
\begin{bmatrix}
\Pr( V_1 , \neg V_2, \neg V_3, S) + \Pr( \neg V_1 , V_2, \neg V_3, S) + \Pr( V_1 , \neg V_2, \neg V_3, S) +\ \\
 \Pr( V_1 , V_2, \neg V_3, \neg S) + \Pr( V_1 , \neg V_2, V_3, \neg S) + \Pr( \neg V_1 , V_2, V_3, \neg S)\\
\end{bmatrix}
}{  
\Omega_{V_1, V_2, V_3, S}
}
\\[10pt]
\Pr(\text{one dose})
&=
\frac{ 
\Pr(V_1, \neg V_2, \neg V_3, \neg S) + \Pr(\neg V_1 , V_2, \neg V_3, \neg S) + \Pr( \neg V_1 , \neg V_2, V_3, \neg S) + \Pr( \neg V_1 , \neg V_2, \neg V_3, S)
}{  
\Omega_{V_1, V_2, V_3, S}
}
\\[10pt]
\Pr(\text{zero doses})
&=
\frac{
\Pr( \neg V_1 , \neg V_2, \neg V_3, \neg S)
}{  
\Omega_{V_1, V_2, V_3, S}
}
\end{aligned}
$$

Using our definition of $\Pr(\text{prior dose})$ for the three-dose model above we can calculate the probability of receiving a dose through an SIA campaign conditionaed on prior vaccination by defining the dropout rate from prior vaccination to SIA campaign: 
$$
d_S = \begin{cases}
          0, & \text{if} \ S \ge \Pr(\text{prior dose}). \\[5pt]
          \frac{(\Pr(\text{prior dose}) - S)}{\Pr(\text{prior dose})}, & \text{if} \ S < \Pr(\text{prior dose}).
          \end{cases}
$$

The expectation for the conditional probability terms in the three-dose model with SIA campaign can then be calculated as shown below. Since there are many terms to define, they are organized by doses received.

Possible ways to receive four doses:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( V_1 , V_2, V_3, S) \big] &= \Pr(V_3 \mid V_1,V_2) \Pr( V_2 \mid V_1) \Pr( V_1 ) \Pr(S \mid \text{prior dose}) \\ &=V_1(1-d_{1,2})(1-d_{2,3})(1-d_S) 
\end{aligned}
$$

Possible ways to receive three doses:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( V_1 , V_2, V_3, \neg S) \big] &= \Pr(V_3 \mid V_1,V_2) \Pr( V_2 \mid V_1) \Pr( V_1 ) \Pr(\neg S \mid \text{prior dose}) \\ &=V_1(1-d_{1,2})(1-d_{2,3})d_S 
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , V_2, \neg V_3, S) \big] &= \Pr(\neg V_3 \mid V_1,V_2) \Pr( V_2 \mid V_1) \Pr( V_1 ) \Pr(S \mid \text{prior dose}) \\ 
&= V_1(1-d_{1,2})d_{2,3}(1-d_S)  
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , \neg V_2, V_3, S) \big] &= \Pr(V_3 \mid V_1, \neg V_2) \Pr(\neg V_2 \mid V_1) \Pr(V_1) \Pr(S \mid \text{prior dose}) \\ 
&= 
\begin{cases}
     0,                                                        & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     V_1 d_{1,2} \big(V_3 - (V_1(1-d_{1,2})) \big)(1-d_S) ,    & \text{if} \ V_3 > \Pr(V_2, V_1).
\end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , V_2, V_3, S) \big] &= \Pr(V_3 \mid \neg V_1, V_2) \Pr( V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr(S \mid \text{prior dose}) \\
&= 
\begin{cases}
     0,                                                               & \text{if} \ V_2 \leq V_1. \\
     0,                                                               & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     \big((1-V_1)(V_2-V_1)\big)  \big(V_3 - (V_1(1-d_{1,2})) \big)(1-d_S),   & \text{otherwise}.
\end{cases}
\\[5pt]
\end{aligned}
$$

Possible ways to receive two doses:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( V_1 , V_2, \neg V_3, \neg S) \big] &= \Pr(\neg V_3 \mid V_1,V_2) \Pr( V_2 \mid V_1) \Pr( V_1 ) \Pr(\neg S \mid \text{prior dose}) \\ 
&= V_1(1-d_{1,2})d_{2,3}d_S
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , \neg V_2, V_3, \neg S) \big] &= \Pr(V_3 \mid V_1, \neg V_2) \Pr(\neg V_2 \mid V_1) \Pr(V_1) \Pr(\neg S \mid \text{prior dose}) \\ 
&= 
\begin{cases}
     0,                                                    & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     V_1 d_{1,2} \big(V_3 - (V_1(1-d_{1,2})) \big)d_S,     & \text{if} \ V_3 > \Pr(V_2, V_1).
\end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , V_2, V_3, \neg S) \big] &= \Pr(V_3 \mid \neg V_1, V_2) \Pr( V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr(\neg S \mid \text{prior dose}) \\
&= 
\begin{cases}
     0,                                                                  & \text{if} \ V_2 \leq V_1. \\
     0,                                                                  & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     \big((1-V_1)(V_2-V_1)\big)  \big(V_3 - (V_1(1-d_{1,2})) \big)d_S,   & \text{otherwise}.
\end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( V_1 , \neg V_2, \neg V_3, S) \big] &= \Pr(\neg V_3 \mid V_1, \neg V_2) \Pr( \neg V_2 \mid V_1) \Pr(V_1) \Pr(S \mid \text{prior dose}) \\
&= 
\begin{cases}
     V_1 d_{1,2}(1-d_S),                                                & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     V_1 d_{1,2} \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big)(1-d_S),        & \text{if} \ V_3 > \Pr(V_2, V_1). 
\end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , V_2, \neg V_3, S) \big] &=
\Pr(\neg V_3 \mid \neg V_1, V_2) \Pr( V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr(S \mid \text{prior dose}) \\
&= 
\begin{cases}
     0,                                                                       & \text{if} \ V_2 \leq V_1. \\
     0,                                                                       & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     (1-V_1)(V_2-V_1) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big) (1-d_S),        & \text{otherwise}.
\end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , \neg V_2, V_3, S) \big] &= \Pr(V_3 \mid \neg V_1, \neg V_2) \Pr( \neg V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr(S \mid \text{prior dose}) \\
&= 
     \begin{cases}
      0,                                                                                                  & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
      (1-V_1) \big( V_3 - V_1 d_{1,2} - (1-V_1)(V_2-V_1) - V_1(1-d_{1,2}) \big) (1-d_S),                           & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\
      (1-V_1) \big(1-(V_2-V_1)\big) \big( V_3 - V_1 d_{1,2} - (1-V_1)(V_2-V_1) - V_1(1-d_{1,2}) \big) (1-d_S),     & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 > V_1.
     \end{cases}
\\[5pt]
\end{aligned}
$$

Possible ways to receive one dose:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( V_1 , \neg V_2, \neg V_3, \neg S) \big] &= \Pr(\neg V_3 \mid V_1, \neg V_2) \Pr( \neg V_2 \mid V_1) \Pr(V_1) \Pr(\neg S \mid \text{prior dose}) \\
&= 
\begin{cases}
     V_1 d_{1,2} d_S,                                                   & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     V_1 d_{1,2} \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big) d_S,           & \text{if} \ V_3 > \Pr(V_2, V_1). 
\end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , V_2, \neg V_3, \neg S) \big] &=
\Pr(\neg V_3 \mid \neg V_1, V_2) \Pr( V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr( \neg S \mid \text{prior dose}) \\
&= 
\begin{cases}
     0,                                                                       & \text{if} \ V_2 \leq V_1. \\
     0,                                                                       & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
     (1-V_1)(V_2-V_1) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big) d_S,            & \text{otherwise}.
\end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , \neg V_2, V_3, \neg S) \big] &= \Pr(V_3 \mid \neg V_1, \neg V_2) \Pr( \neg V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr( \neg S \mid \text{prior dose}) \\
&= 
     \begin{cases}
      0,                                                                                                  & \text{if} \ V_3 \leq \Pr(V_2, V_1). \\
      (1-V_1) \big( V_3 - V_1 d_{1,2} - (1-V_1)(V_2-V_1) - V_1(1-d_{1,2}) \big) d_S,                           & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\
      (1-V_1) \big(1-(V_2-V_1)\big) \big( V_3 - V_1 d_{1,2} - (1-V_1)(V_2-V_1) - V_1(1-d_{1,2}) \big) d_S,     & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 > V_1.
     \end{cases}
\\[5pt]
\mathbb{E}\big[ \Pr( \neg V_1 , \neg V_2, \neg V_3, S) \big] &= \Pr(\neg V_3 \mid \neg V_1, \neg V_2) \Pr( \neg V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr(S \mid \text{prior dose})  \\
&= 
\begin{cases}
      (1-V_1)(1-d_S),                                                                    & \text{if} \ V_3 \leq \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\
      (1-V_1) \big(1-(V_2-V_1)\big)(1-d_S),                                            & \text{if} \ V_3 \leq \Pr(V_2, V_1) \ \& \ V_2 >  V_1. \\
      (1-V_1) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big)(1-d_S),                          & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\
      (1-V_1) \big(1-(V_2-V_1)\big) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big)(1-d_S),    & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 >  V_1.
\end{cases}
\\[5pt]
\end{aligned}
$$

Possible ways to receive zero doses:
$$
\begin{aligned}
\mathbb{E}\big[ \Pr( \neg V_1 , \neg V_2, \neg V_3, \neg S) \big] &= \Pr(\neg V_3 \mid \neg V_1, \neg V_2) \Pr( \neg V_2 \mid \neg V_1) \Pr(\neg V_1) \Pr(\neg S \mid \text{prior dose})  \\
&= 
\begin{cases}
      (1-V_1)d_S,                                                                    & \text{if} \ V_3 \leq \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\
      (1-V_1) \big(1-(V_2-V_1)\big)d_S,                                            & \text{if} \ V_3 \leq \Pr(V_2, V_1) \ \& \ V_2 >  V_1. \\
      (1-V_1) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big)d_S,                          & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 \leq V_1. \\
      (1-V_1) \big(1-(V_2-V_1)\big) \big(1 - (V_3 - (V_1(1-d_{1,2}))) \big)d_S,    & \text{if} \ V_3 >  \Pr(V_2, V_1) \ \& \ V_2 >  V_1.
\end{cases}
\end{aligned}
$$


### Total proportion immune due to vaccination
Given the definitions of receiving 0, 1, 2, 3, and 4 doses shown above, the total proportion vaccinated $p_\text{vacc}$ under the three-dose model with an SIA:
$$
p_{\text{vacc}} = \sum_j \text{effectiveness}_j \times \Pr(\text{dose}_j), \ \  j = \{ 1, 2, 3, 4\}
$$


### Caveats

  * Modeling the proportion of the population that has received $k$ doses of a vaccine can help to make better estimates of the amount of population immunity due to vaccination. However, when infection imbues temporary or life-long immunity, the total proportion immune will include those that have been exposed or infected.
Therefore to more fully estimate the total proportion immune $p_{\text{imm}}$ we would have to estimate the sum of the proportion vaccinated and the portion infected and/or recovered.
$$
p_{\text{imm}} = p_{\text{vacc}} + p_{\text{inf}}
$$

  * The method does not explicitly account for age groups, but it can be applied to each age group for which vaccination coverage data are available.
  
  * Does not explicitly include vaccine effectiveness or efficacy, however effectiveness can be included in the definition of $V_k$, and efficacy can be included in the estimation of total proportion vaccinated.
  
  * Models include 2-dose routine, 2-dose routine + SIA, 3-dose routine, and 3-dose routine + SIA, however, vaccination schedules vary across countries. Formulations that are easily expandable to larger numbers of routine and campaign activities assume independence among activities. More difficult to generalize condition probabilities to any number of vaccination events.
  
  * Previous methods that assume independence among routine activities and SIAs likely over-estimate immunity. The dependence assumptions here are conservative, so the estimates of immunity will tend to be under-estimated.


