23a24
>         real8   timetostop;
49a51
>         
54a57
>         timetostop		 = home->param->stopTime;
55a59,60
>         
>         
93c98
< 
---
> param->velstatTime=param->timeNow;
100,101c105,108
< 
<         while (home->cycle < cycleEnd) {
---
> /* (home->cycle < cycleEnd) &&*/ 
>         while ((home->cycle < cycleEnd) && (param->timeNow<timetostop)) {
> 			/*printf("%e %e", param->timeNow,timetostop);*/
> 			
102a110
>            
