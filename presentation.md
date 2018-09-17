# Presentation

## Prologue
### Greeting
Hello guys, thank you for comming. Today I'm gonna show you an awesome topic regarding how to price and hedge American-style options. Really a well-defined and classic problem, right? So yeah, with all of the new and improved techniques, that is our project.

### Gratitude
At the start of my report, I want to first appreciate Prof. Cui for his continuous guidance and interest on my research. And also I want to thank Wan, Tengxi, one of my friend, some of you may know him as he took this course last semester and my research project is kind of extension based on his interests on LSM.

## Introduction
### Outline
OK, back to our report. Today, I will go through some basic theoretical details about the objective problem, and then, show you some stunning methods to address it. In addition, select several potential extentions or research directions to further hone this project.

### Background
As we know, option contracts play a prominent role in financial industry including equity, fixed income, forex, commodities and economic index. It can be divided into two main types by exercise rule: European options and American options. Today, We are going to focus on American options and how to price and hedge it step by step. I know the subject today is literature review, but I don't want to be boring, how can anybody convey the ideas without a proper clarification.

### Pricing
Here is the pricing formula under the optimization problem. For the sake of simplicity, I set the risk free rate as 0 here, otherwise, one can view it as a forward price standing at the inception. Note that the price process is a supermartingale. (intuitively, the earlier you are, the more chance to find a better exercise point)

### Hedging
Besides pricing, hedging is increasingly important in nowadays especially for the ubiquitous American options. But the question is there is no explicite formula for such a option. That means, we cannot obtain the hedging strategies by directly differential. Here, Doob-Meyer decompostion of the supermartingle representation stands out and find a way to calculate hedging strategy corresponding to the underlying asset.

### Phase transition
In particular, there comes a question, if we assert the American options' price is actually a supermartingale instead of a martingale, is it contrary to the martingale pricing theory? Well, from I'm concerned, it's not. Before the underlying arrive the exercise boundary($s\in[t,\tau]$),  the price process is a martingale. But after that, it's a supermartingale as you lose your chance in perspective of probability (you're not supposed to do so). Fair enough, right? Correct me if I'm wrong.

### Methods Classification
In this project, we focuse on Monte Carlo simulation method, as it has many advantages comparing to others I will tell you in the following parts. In a recent research conducted by Bouchard and Warin, they divid Monte Carlo simulation method in the problem of pricing and hedging into two approaches: Malliavin-based approach and Regression-based approach.

### Purpose
As I mentioned, this project focuses on the Monte Carlo method. In order to make it more original and profound, it's not enough to just replicate people's work, there are several ideas I want to tackle. Basically, just a rough idea, and I wish in a foreseeable future we can achieve them.

## Regression Based Methods

### LSM
* Longstaff & Schiwartz (2001)

### LSM improvements
* CIR (1985)
* Anderson (2007)
* Heston (1993)

### OHMC
* Potters, Bouchaud & Sestovic (2001)

### Other extensions
* Bouchard & Warin (2012)
* Davis & Zariphopoulou (1995)
* Broadie(2011)

## Conclusion
pros:

* Simple(algorithm) yet powerfull（expansibility)
* Compatible to any model (nonparametric)
* curse of dimension

cons:

* bias
* computational complexity

## Thank you

These are my references upon now, in case you want to follow this topic. So, yeah, that's it. Thank you for attention. Any question?

