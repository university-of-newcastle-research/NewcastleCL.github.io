---
layout: post
title: PMWG version 0.2.0 released
---

The second CRAN release of the PMWG package is now available for download and includes some new functionality, improvements in the sampling stage and some minor bugfixes. I'll cover the major points in this post, but further details are available in the changelog for the package (https://newcastlecl.github.io/pmwg/news/index.html).

## Main improvements

### Updating the proposal distribution

One of the biggest parts of this updated version is the addition of a new process for updating the proposal distribution during the sampling stage. The proposal distribution previously was conditioned on the samples from the adaptation stage and led to major efficiency gains in the final sampling stage. The latest version will now update this conditional distribution using samples from both adaptation and the sampling stage every 500 iterations.

The primary result of this is that if a proposal distribution at the end of the adaptation stage is not using samples from an appropriate region of the posterior then it will not be as efficient as possible. Therefore updating it during the sampling stage means there will be a higher probability that the recomputed proposal distribution is conditioned on appropriate samples from the posterior, and will be more efficient.

### Other Changes

There is now a new function, `relabel_samples` that can change the designation of the source stage for samples in the pmwgs object. The purpose of this is to enable the end user to reallocate burn in samples to adaptation samples in the case where sampling is expensive and a subset of the first stage samples are sufficiently burnt in. The samples can be relabelled as adaptation samples using this function, and can then be used to inform the proposal distribution.

The display of the acceptance rate (now labelled **New**) has been updated to show the runnning average of newly accepted particles over the last 200 samples. This means the acceptance rate should be more responsive to the current state of the sampling efficiency.

Finally, there is a new object included as a data structure with the package, called `sampled_forstmann`. This will allow people unfamiliar with the project to trial running more of the API before adopting the package for their own research, as well as allowing better internal testing of the API going forward.

## Other news

In other news it has been exciting to see the interest in the community when presenting talks about the package, and the feedback and questions on Twitter and via email.

A bug thanks to Russell J. Boag from the *Integrative Model-based Cognitive Neuroscience Research Unit* at the University of Amsterdam for his blog post about fitting joint models with the PMwG sampler, which can be found [here](https://github.com/russell-j-boag/russell-j-boag.github.io/blob/main/tutorial_joint_ddm_pmwg.md)


