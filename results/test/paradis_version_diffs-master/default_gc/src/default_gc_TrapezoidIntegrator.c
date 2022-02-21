8a9
> #include "Comm.h"
596a598
>                
745c747
<                 if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
---
>                 if ((newDT < 1.0e-25) && (home->myDomain == 0)) {
