#!/bin/bash
# NOTE (after changing shell flag in contact_run.R to TRUE)
# arg1: sensitivity analysis: none, base/lo/hi tptru (SOC resource use), hicoprev (higher co-prevalence based on systematic review)

R --slave --vanilla --args <contact_run.R lo & R --slave --vanilla --args <contact_run.R hi & R --slave --vanilla --args <contact_run.R tptru & R --slave --vanilla --args <contact_run.R hicoprev
R --slave --vanilla --args <contact_run.R none



