 /*-- renderGraph $Revision: 1.3 $ -- 
 -- parameters $Revision: 1.5 $ -- 
 -- parameterUse $Revision: 1.2 $ -- 
 -- arrange $Revision: 1.2 $ -- 
 -- arrangePlus $Revision: 1.2 $ -- 
 -- arrangeAnnotate $Revision: 1.2 $ -- 
// generating...hept
 -- arrangeMain $Revision: 1.3 $ -- 
 -- arrangeInvariant $Revision: 1.2 $ -- 
// stack sizes = 13 3
// cases considered= 1000
// stack sizes = 11 6
// cases considered= 2000
// stack sizes = 12 8
// cases considered= 3000
// stack sizes = 13 10
// cases considered= 4000
// stack sizes = 7 10
// cases considered= 5000
// stack sizes = 6 16
// cases considered= 6000
// stack sizes = 16 16
// cases considered= 7000
// stack sizes = 12 16
// cases considered= 8000
// stack sizes = 20 16
// cases considered= 9000
// stack sizes = 8 18
// cases considered= 10000
// stack sizes = 15 20
// cases considered= 11000
// stack sizes = 19 20
// cases considered= 12000
// stack sizes = 17 20
// cases considered= 13000
// stack sizes = 9 23
// cases considered= 14000
// stack sizes = 7 25
// cases considered= 15000
// stack sizes = 27 25
// cases considered= 16000
// stack sizes = 21 25
// cases considered= 17000
// stack sizes = 15 25
// cases considered= 18000
// stack sizes = 20 25
// cases considered= 19000
// stack sizes = 16 25
// cases considered= 20000
// stack sizes = 5 25
// cases considered= 21000
// stack sizes = 19 26
// cases considered= 22000
// stack sizes = 6 26
// cases considered= 23000
// stack sizes = 8 26
// cases considered= 24000
// stack sizes = 23 26
// cases considered= 25000
// stack sizes = 15 26
// cases considered= 26000
// stack sizes = 3 26
// cases considered= 27000
// stack sizes = 5 26
// cases considered= 28000
// stack sizes = 22 26
// cases considered= 29000
// stack sizes = 7 26
// cases considered= 30000
// stack sizes = 5 26
// cases considered= 31000
// Total pops: 31885
stack size = 26
count = 31885
*/
public class graphHEPT {
final static String data[] = {

" 0 15  7 8 11 7 1 0 2 3 4 8 3 0 4 3 3 2 0 3 4 0 5 3 5 0 1 3 5 1 6 3 6 1 7 3 6 7 10 3 10 7 11 3 10 11 8 3 10 8 9 3 9 8 4 3 6 10 9 3 5 6 9 3 9 4 5",

" 0 15  7 9 11 7 2 0 3 8 4 9 8 3 4 3 4 3 0 3 4 0 1 3 1 0 2 3 1 2 6 3 6 2 7 3 6 7 10 3 10 7 11 3 10 11 9 3 10 9 5 3 5 9 4 3 5 4 1 3 5 1 6 3 6 10 5",

" 0 15  7 8 4 0 1 5 9 11 4 8 11 9 10 3 10 9 5 3 10 5 6 3 6 5 2 3 2 5 1 3 2 1 0 3 2 0 3 3 6 2 3 3 3 0 4 3 3 4 7 3 6 3 7 3 10 6 7 3 8 10 7 3 7 4 8",

" 0 15  7 8 5 1 2 6 9 11 4 8 11 10 7 3 10 11 9 3 10 9 6 3 7 10 6 3 7 6 3 3 3 6 2 3 3 2 0 3 0 2 1 3 0 1 5 3 0 5 4 3 3 0 4 3 7 3 4 3 8 7 4 3 4 5 8",

" 0 15  7 0 1 5 10 11 9 4 4 0 4 8 3 3 8 4 9 3 8 9 11 3 8 11 7 3 3 8 7 3 7 11 10 3 7 10 6 3 3 7 6 3 6 10 5 3 6 5 2 3 3 6 2 3 0 3 2 3 1 0 2 3 2 5 1",

" 0 15  7 7 10 11 9 6 2 3 4 7 3 0 4 3 0 3 2 3 0 2 1 3 1 2 6 3 1 6 5 3 5 6 9 3 5 9 8 3 8 9 11 3 8 11 7 3 11 10 7 3 8 7 4 3 5 8 4 3 1 5 4 3 4 0 1",

" 0 17  7 2 9 11 6 7 12 8 3 2 8 7 3 8 12 7 4 2 7 1 0 3 1 7 6 3 1 6 5 3 5 6 11 3 5 11 10 3 10 11 9 3 10 9 3 3 3 9 2 3 3 2 0 3 3 0 4 3 10 3 4 3 5 10 4 3 1 5 4 3 4 0 1",

" 0 17  7 0 1 2 3 4 5 6 3 0 6 5 4 0 5 4 7 3 7 4 8 3 8 4 3 3 8 3 9 3 9 3 2 3 9 2 10 3 10 2 1 3 10 1 11 3 11 1 0 3 11 0 7 3 11 7 12 3 12 7 8 3 12 8 9 3 12 9 10 3 10 11 12",

" 0 17  7 9 12 7 2 0 3 8 3 9 8 3 4 9 3 4 10 3 4 3 0 3 4 0 1 3 1 0 2 3 1 2 6 3 6 2 7 3 6 7 11 3 11 7 12 3 11 12 10 3 12 9 10 3 11 10 5 3 5 10 4 3 5 4 1 3 5 1 6 3 6 11 5",

" 0 15  7 6 10 11 9 4 1 5 3 6 5 1 4 6 1 0 2 3 0 1 4 3 0 4 3 3 3 4 9 3 3 9 8 3 8 9 11 3 8 11 7 3 7 11 10 3 7 10 6 3 7 6 2 3 8 7 2 3 3 8 2 3 2 0 3",

" 0 14  7 8 5 1 2 6 9 10 3 8 10 9 3 8 9 7 3 7 9 6 3 7 6 3 3 3 6 2 3 3 2 0 3 0 2 1 3 0 1 5 3 0 5 4 3 4 5 8 3 4 8 7 3 3 0 4 3 4 7 3",

" 0 16  7 7 8 9 10 5 6 11 3 7 11 6 3 7 6 1 3 1 6 5 3 1 5 0 3 0 5 4 3 4 5 10 3 4 10 9 3 4 9 3 3 3 9 8 3 3 8 2 3 2 8 7 3 2 7 1 3 0 4 3 3 0 3 2 3 2 1 0",

" 0 18  7 3 0 4 11 12 9 10 3 3 10 9 3 3 9 8 3 8 9 12 3 8 12 7 3 7 12 6 3 6 12 11 3 6 11 5 3 5 11 4 3 5 4 0 3 5 0 1 3 6 5 1 3 1 0 2 3 2 0 3 3 2 3 8 3 2 8 7 3 7 6 1 3 1 2 7",

" 0 14  7 1 3 8 10 7 2 0 3 1 0 2 3 1 2 5 3 5 2 6 3 6 2 7 3 6 7 10 3 6 10 9 3 9 10 8 3 9 8 4 3 4 8 3 3 4 3 1 3 4 1 5 3 5 6 9 3 9 4 5",

" 0 16  7 1 6 10 9 4 0 5 3 1 5 0 3 1 0 2 3 2 0 3 3 3 0 4 3 3 4 9 3 3 9 8 3 8 9 11 3 11 9 10 3 11 10 6 3 11 6 7 3 8 11 7 3 7 6 1 3 7 1 2 3 2 3 8 3 8 7 2",

" 0 18  7 1 0 2 8 9 6 7 3 1 7 6 3 1 6 5 3 5 6 12 3 12 6 9 3 12 9 10 3 10 9 2 3 9 8 2 3 10 2 3 3 3 2 0 3 3 0 4 3 4 0 1 3 4 1 5 3 4 5 11 3 11 5 12 3 11 12 10 3 11 10 3 3 3 4 11",

" 0 15  7 8 5 1 2 6 9 11 3 8 11 10 3 10 11 9 3 10 9 6 4 8 10 6 7 3 7 6 3 3 3 6 2 3 3 2 0 3 0 2 1 3 0 1 5 3 0 5 4 3 3 0 4 3 7 3 4 3 8 7 4 3 4 5 8",

" 0 19  7 0 1 2 3 4 5 6 3 0 6 7 3 7 6 5 3 7 5 4 4 0 7 4 8 3 8 4 9 3 9 4 3 3 9 3 10 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 12 3 12 1 0 3 12 0 8 3 12 8 13 3 13 8 9 3 13 9 10 3 13 10 11 3 11 12 13",

" 0 14  7 10 8 3 0 1 4 9 3 10 9 5 3 5 9 4 3 5 4 1 3 5 1 6 3 6 1 2 3 2 1 0 3 2 0 3 3 2 3 7 3 7 3 8 3 7 8 10 3 7 10 6 3 10 5 6 3 6 2 7",

" 0 16  7 6 1 2 3 7 11 10 3 6 10 9 3 9 10 11 3 9 11 7 3 9 7 8 3 8 7 4 3 4 7 3 3 4 3 0 3 0 3 2 3 0 2 1 3 0 1 5 3 4 0 5 3 5 1 6 3 5 6 8 3 6 9 8 3 8 4 5",

" 0 15  7 10 4 5 6 7 8 11 3 10 11 9 3 9 11 8 3 9 8 2 3 2 8 7 3 2 7 1 3 1 7 6 3 1 6 5 3 1 5 0 3 0 5 4 3 0 4 3 3 3 4 10 3 3 10 9 4 0 3 9 2 3 2 1 0",

" 0 15  7 0 1 2 3 4 5 6 3 0 6 7 3 7 6 5 3 7 5 8 3 8 5 4 3 8 4 9 3 9 4 3 3 9 3 10 3 10 3 2 3 10 2 1 3 10 1 0 4 9 10 0 11 3 11 0 7 3 11 7 8 3 8 9 11",

" 0 15  7 7 11 10 6 1 0 2 3 7 2 3 3 3 2 0 3 3 0 4 3 4 0 1 3 4 1 5 3 5 1 6 3 5 6 9 3 9 6 10 3 9 10 11 3 9 11 8 3 8 11 7 4 8 7 3 4 3 8 4 5 3 5 9 8",

" 0 17  7 0 1 2 3 4 5 6 3 0 6 7 3 7 6 5 3 7 5 8 3 8 5 4 3 8 4 9 3 9 4 3 3 9 3 10 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 4 11 0 7 12 3 12 7 8 3 12 8 9 3 12 9 10 3 10 11 12",

" 0 17  7 0 1 2 3 4 5 6 3 0 6 7 3 7 6 5 3 7 5 8 3 8 5 4 3 8 4 9 3 9 4 3 3 9 3 10 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 4 10 11 0 12 3 12 0 7 3 12 7 8 3 12 8 9 3 9 10 12",

" 0 19  7 7 12 11 6 2 0 3 3 7 3 4 3 4 3 1 3 1 3 0 3 1 0 2 3 1 2 5 3 5 2 6 3 5 6 10 3 10 6 11 3 10 11 13 3 13 11 12 3 13 12 8 3 8 12 7 3 8 7 4 3 8 4 9 4 9 4 1 5 3 9 5 10 3 9 10 13 3 13 8 9"};

};
// Generator complete. All 7-gons determined
