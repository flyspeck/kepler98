 /*-- renderGraph $Revision: 1.3 $ -- 
 -- parameters $Revision: 1.5 $ -- 
 -- parameterUse $Revision: 1.2 $ -- 
 -- arrange $Revision: 1.2 $ -- 
 -- arrangePlus $Revision: 1.2 $ -- 
 -- arrangeAnnotate $Revision: 1.2 $ -- 
// generating...oct
 -- arrangeMain $Revision: 1.3 $ -- 
 -- arrangeInvariant $Revision: 1.2 $ -- 
// stack sizes = 13 1
// cases considered= 1000
// stack sizes = 11 1
// cases considered= 2000
// stack sizes = 8 6
// cases considered= 3000
// stack sizes = 11 6
// cases considered= 4000
// stack sizes = 16 9
// cases considered= 5000
// stack sizes = 8 14
// cases considered= 6000
// stack sizes = 5 16
// cases considered= 7000
// stack sizes = 2 18
// cases considered= 8000
// stack sizes = 18 18
// cases considered= 9000
// stack sizes = 12 18
// cases considered= 10000
// stack sizes = 9 19
// cases considered= 11000
// stack sizes = 23 22
// cases considered= 12000
// stack sizes = 18 22
// cases considered= 13000
// stack sizes = 4 25
// cases considered= 14000
// stack sizes = 5 26
// cases considered= 15000
// stack sizes = 13 28
// cases considered= 16000
// stack sizes = 16 28
// cases considered= 17000
// stack sizes = 12 28
// cases considered= 18000
// stack sizes = 5 28
// cases considered= 19000
// stack sizes = 13 28
// cases considered= 20000
// stack sizes = 7 28
// cases considered= 21000
// stack sizes = 18 28
// cases considered= 22000
// stack sizes = 12 28
// cases considered= 23000
// stack sizes = 25 29
// cases considered= 24000
// stack sizes = 5 29
// cases considered= 25000
// stack sizes = 10 33
// cases considered= 26000
// stack sizes = 5 33
// cases considered= 27000
// stack sizes = 8 33
// cases considered= 28000
// stack sizes = 10 34
// cases considered= 29000
// stack sizes = 6 34
// cases considered= 30000
// stack sizes = 17 34
// cases considered= 31000
// stack sizes = 5 34
// cases considered= 32000
// stack sizes = 5 34
// cases considered= 33000
// stack sizes = 6 34
// cases considered= 34000
// stack sizes = 8 34
// cases considered= 35000
// stack sizes = 7 34
// cases considered= 36000
// stack sizes = 15 34
// cases considered= 37000
// stack sizes = 7 34
// cases considered= 38000
// stack sizes = 8 34
// cases considered= 39000
// stack sizes = 6 34
// cases considered= 40000
// stack sizes = 1 34
// cases considered= 41000
// Total pops: 41573
stack size = 34
count = 41573
*/
public class graphOCT {
final static String data[] = {

" 0 14  8 0 1 2 3 4 5 6 7 4 0 7 6 8 3 8 6 9 3 9 6 5 3 9 5 4 3 9 4 10 3 8 9 10 3 10 4 3 3 10 3 11 3 8 10 11 3 0 8 11 3 11 3 2 3 0 11 2 3 2 1 0",

" 0 13  8 8 6 1 2 3 7 9 10 3 8 10 7 3 10 9 7 3 8 7 4 3 4 7 3 3 4 3 0 3 0 3 2 3 0 2 1 3 0 1 5 3 4 0 5 3 5 1 6 3 5 6 8 3 8 4 5",

" 0 15  8 8 5 1 2 6 9 11 10 3 8 10 9 3 10 11 9 3 8 9 7 3 7 9 6 3 7 6 3 3 3 6 2 3 3 2 0 3 0 2 1 3 0 1 5 3 0 5 4 3 3 0 4 3 7 3 4 3 8 7 4 3 4 5 8",

" 0 17  8 1 4 8 12 7 2 3 0 3 1 0 2 3 0 3 2 3 1 2 6 3 6 2 7 3 6 7 10 3 10 7 11 3 11 7 12 3 11 12 8 3 11 8 9 3 10 11 9 3 9 8 4 3 9 4 5 3 10 9 5 3 6 10 5 3 1 6 5 3 5 4 1",

" 0 19  8 8 3 9 10 6 11 13 12 3 8 12 11 3 12 13 11 3 8 11 7 3 7 11 6 3 7 6 1 3 1 6 5 3 5 6 10 3 5 10 4 3 4 10 9 3 4 9 3 3 4 3 0 3 5 4 0 3 1 5 0 3 0 3 2 3 2 3 8 3 2 8 7 3 2 7 1 3 1 0 2",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 6 4 0 6 5 8 3 8 5 9 3 9 5 4 3 9 4 3 3 9 3 10 3 8 9 10 3 10 3 2 3 10 2 11 3 8 10 11 3 0 8 11 3 1 0 11 3 11 2 1",

" 0 13  8 4 8 10 7 2 0 1 3 3 4 3 1 3 4 1 5 3 5 1 2 3 1 0 2 3 5 2 6 3 6 2 7 3 6 7 10 3 6 10 9 3 5 6 9 3 9 10 8 3 9 8 4 3 4 5 9",

" 0 15  8 6 1 2 3 7 11 9 10 3 6 10 9 3 6 9 8 3 8 9 7 3 9 11 7 3 8 7 4 3 4 7 3 3 4 3 0 3 0 3 2 3 0 2 1 3 0 1 5 3 4 0 5 3 8 4 5 3 6 8 5 3 5 1 6",

" 0 13  8 7 8 9 4 0 5 6 10 3 7 10 6 3 7 6 1 3 1 6 5 3 1 5 0 3 1 0 2 3 2 0 3 3 3 0 4 3 3 4 9 3 3 9 8 3 2 3 8 3 7 1 2 3 2 8 7",

" 0 13  8 6 10 8 3 4 0 5 9 3 6 9 5 3 6 5 1 3 1 5 0 3 1 0 2 3 2 0 3 3 0 4 3 3 2 3 8 3 2 8 7 3 7 8 10 3 7 10 6 3 7 6 1 3 1 2 7",

" 0 15  8 2 5 10 8 9 4 1 0 3 2 0 1 3 2 1 3 3 3 1 4 3 3 4 7 3 7 4 8 3 4 9 8 3 7 8 11 3 11 8 10 3 11 10 5 3 11 5 6 3 7 11 6 3 6 5 2 3 6 2 3 3 3 7 6",

" 0 13  8 4 8 9 10 7 2 0 3 3 4 3 0 3 4 0 1 3 1 0 2 3 1 2 6 3 6 2 7 3 6 7 9 3 7 10 9 3 6 9 5 3 5 9 8 3 5 8 4 3 5 4 1 3 1 6 5",

" 0 15  8 0 4 8 11 10 7 2 3 3 0 3 2 3 0 2 1 3 1 2 7 3 1 7 6 3 6 7 10 3 6 10 9 3 9 10 8 3 10 11 8 3 9 8 4 3 9 4 5 3 5 4 0 3 5 0 1 3 5 1 6 3 6 9 5",

" 0 15  8 0 1 6 7 11 10 4 5 3 0 5 4 3 0 4 3 3 3 4 10 3 3 10 9 3 9 10 11 3 9 11 8 3 8 11 7 3 8 7 1 3 7 6 1 3 8 1 2 3 2 1 0 3 2 0 3 3 2 3 9 3 9 8 2",

" 0 17  8 1 3 7 12 11 6 2 0 3 1 0 2 3 1 2 5 3 5 2 6 3 5 6 10 3 10 6 11 3 10 11 9 3 9 11 8 3 8 11 12 3 8 12 7 3 8 7 3 3 9 8 3 3 9 3 4 3 4 3 1 3 4 1 5 3 4 5 10 3 10 9 4",

" 0 19  8 10 6 11 12 8 3 9 13 3 10 13 9 3 10 9 4 3 4 9 3 3 4 3 0 3 0 3 2 3 2 3 8 3 2 8 7 3 7 8 11 3 8 12 11 3 7 11 6 3 7 6 1 3 2 7 1 3 1 6 5 3 5 6 10 3 5 10 4 3 5 4 0 3 0 2 1 3 1 5 0",

" 0 13  8 4 1 0 2 5 10 8 9 3 4 9 8 3 4 8 7 3 7 8 6 3 6 8 10 3 6 10 5 3 6 5 2 3 7 6 2 3 7 2 3 3 3 2 0 3 3 0 1 3 3 1 4 3 4 7 3",

" 0 15  8 10 7 2 0 1 3 8 11 3 10 11 8 3 10 8 9 3 9 8 4 3 4 8 3 3 4 3 1 3 4 1 5 3 5 1 2 3 1 0 2 3 5 2 6 3 6 2 7 3 9 4 5 3 9 5 6 3 10 9 6 3 6 7 10",

" 0 13  8 3 6 10 9 5 2 0 1 3 3 1 4 3 4 1 2 3 1 0 2 3 4 2 5 3 4 5 8 3 8 5 9 3 8 9 7 3 7 9 10 3 7 10 6 3 7 6 3 3 8 7 3 3 3 4 8",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 5 3 7 6 5 3 8 5 4 3 8 4 9 3 9 4 3 3 9 3 10 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 4 9 10 11 0 3 0 8 9",

" 0 15  8 5 9 11 8 3 0 4 1 3 5 1 2 3 2 1 0 3 1 4 0 3 2 0 3 3 2 3 6 3 6 3 7 3 7 3 8 3 7 8 11 3 7 11 10 3 10 11 9 3 10 9 5 3 10 5 6 3 5 2 6 3 6 7 10",

" 0 17  8 8 4 5 1 6 9 12 11 3 8 11 10 3 10 11 9 3 11 12 9 3 10 9 6 3 10 6 7 3 7 6 2 3 2 6 1 3 2 1 0 3 0 1 5 3 0 5 4 3 0 4 3 3 2 0 3 3 8 10 7 3 7 2 3 3 8 7 3 3 3 4 8",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 5 3 8 5 9 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 4 9 10 11 0 3 0 8 9",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 5 3 8 5 9 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 11 3 11 3 2 3 11 2 1 3 11 1 0 4 9 10 11 0 3 0 8 9",

" 0 15  8 8 4 1 0 2 5 9 11 3 8 11 10 3 10 11 9 3 10 9 5 3 10 5 6 3 6 5 2 3 6 2 7 3 7 2 3 3 3 2 0 3 8 10 6 3 3 0 1 3 3 1 4 3 7 3 4 3 8 6 7 3 7 4 8",

" 0 19  8 1 0 2 9 10 12 13 8 3 1 8 7 3 7 8 13 3 7 13 12 3 7 12 6 3 6 12 11 3 11 12 10 3 11 10 4 3 4 10 3 3 3 10 9 3 1 7 6 3 3 9 2 3 3 2 0 3 4 3 0 3 4 0 5 3 5 0 1 3 5 1 6 3 5 6 11 3 11 4 5",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 9 3 9 6 5 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 3 11 0 8 4 8 9 10 11",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 9 3 9 6 5 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 11 3 11 3 2 3 11 2 1 3 11 1 0 4 10 11 0 8 3 8 9 10",

" 0 17  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 5 3 7 6 5 3 9 5 4 3 9 4 10 3 10 4 11 3 11 4 3 3 11 3 2 3 11 2 12 3 10 11 12 3 8 9 10 3 8 10 12 3 0 8 12 3 0 12 2 3 2 1 0",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 4 11 3 0 8 3 0 3 2 3 8 9 10 3 10 11 8 3 2 1 0",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 3 11 3 1 3 3 2 1 3 11 1 0 4 10 11 0 8 3 8 9 10",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 3 11 3 2 4 11 2 0 8 3 2 1 0 3 8 9 10 3 10 11 8",

" 0 16  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 3 11 3 2 3 11 2 12 3 12 2 1 3 12 1 0 3 12 0 8 4 10 11 12 8 3 8 9 10",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 10 3 10 6 5 3 10 5 4 3 10 4 11 3 11 4 3 4 11 3 0 8 3 0 3 2 3 9 10 11 3 11 8 9 3 2 1 0"};

};
// Generator complete. All 8-gons determined
