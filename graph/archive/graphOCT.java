package archive;
/*
Total pops: 54568
stack size = 17
count = 54568
*/
public class graphOCT {
public final static String data[] = {
 
" 0 15  8 11 10 8 5 1 2 6 9 3 10 11 9 3 8 10 9 3 8 9 7 3 7 9 6 3 7 6 3 3 3 6 2 3 3 2 0 3 0 2 1 3 0 1 5 3 0 5 4 3 3 0 4 3 7 3 4 3 8 7 4 3 4 5 8",
 
" 0 17  8 3 0 1 4 8 12 7 2 3 0 3 2 3 1 0 2 3 1 2 6 3 6 2 7 3 6 7 10 3 10 7 11 3 11 7 12 3 11 12 8 3 11 8 9 3 10 11 9 3 9 8 4 3 9 4 5 3 10 9 5 3 6 10 5 3 1 6 5 3 5 4 1",
 
" 0 19  8 13 12 8 3 9 10 6 11 3 12 13 11 3 8 12 11 3 8 11 7 3 7 11 6 3 7 6 1 3 1 6 5 3 5 6 10 3 5 10 4 3 4 10 9 3 4 9 3 3 4 3 0 3 5 4 0 3 1 5 0 3 0 3 2 3 2 3 8 3 2 8 7 3 2 7 1 3 1 0 2",
 
" 0 15  8 10 6 1 2 3 7 11 9 3 6 10 9 3 6 9 8 3 8 9 7 3 9 11 7 3 8 7 4 3 4 7 3 3 4 3 0 3 0 3 2 3 0 2 1 3 0 1 5 3 4 0 5 3 8 4 5 3 6 8 5 3 5 1 6",
 
" 0 15  8 0 2 5 10 8 9 4 1 3 2 0 1 3 2 1 3 3 3 1 4 3 3 4 7 3 7 4 8 3 4 9 8 3 7 8 11 3 11 8 10 3 11 10 5 3 11 5 6 3 7 11 6 3 3 7 6 3 2 3 6 3 6 5 2",
 
" 0 15  8 3 0 4 8 11 10 7 2 3 0 3 2 3 0 2 1 3 1 2 7 3 1 7 6 3 6 7 10 3 6 10 9 3 9 10 8 3 10 11 8 3 9 8 4 3 9 4 5 3 6 9 5 3 1 6 5 3 0 1 5 3 5 4 0",
 
" 0 15  8 5 0 1 6 7 11 10 4 3 0 5 4 3 0 4 3 3 3 4 10 3 3 10 9 3 9 10 11 3 9 11 8 3 8 11 7 3 8 7 1 3 7 6 1 3 8 1 2 3 2 1 0 3 2 0 3 3 2 3 9 3 9 8 2",
 
" 0 17  8 0 1 3 7 12 11 6 2 3 1 0 2 3 1 2 5 3 5 2 6 3 5 6 10 3 10 6 11 3 10 11 9 3 9 11 8 3 8 11 12 3 8 12 7 3 8 7 3 3 9 8 3 3 9 3 4 3 4 3 1 3 4 1 5 3 4 5 10 3 10 9 4",
 
" 0 17  8 11 7 8 9 12 10 5 6 3 7 11 6 3 7 6 1 3 1 6 5 3 1 5 0 3 0 5 4 3 4 5 10 3 4 10 9 3 10 12 9 3 4 9 3 3 3 9 8 3 0 4 3 3 3 8 2 3 0 3 2 3 1 0 2 3 7 1 2 3 2 8 7",
 
" 0 15  8 0 1 3 7 11 10 6 2 3 1 0 2 3 1 2 5 3 5 2 6 3 5 6 8 3 8 6 9 3 9 6 10 3 9 10 11 3 9 11 7 3 8 9 7 3 8 7 4 3 4 7 3 3 5 8 4 3 1 5 4 3 4 3 1",
 
" 0 17  8 11 7 8 12 9 10 5 6 3 7 11 6 3 7 6 1 3 1 6 5 3 1 5 0 3 0 5 4 3 4 5 10 3 4 10 9 3 4 9 3 3 0 4 3 3 3 9 8 3 9 12 8 3 3 8 2 3 2 8 7 3 2 7 1 3 2 1 0 3 0 3 2",
 
" 0 15  8 4 1 5 9 11 8 3 0 3 1 4 0 3 1 0 2 3 2 0 3 3 2 3 6 3 6 3 7 3 7 3 8 3 7 8 11 3 7 11 10 3 6 7 10 3 10 11 9 3 10 9 5 3 6 10 5 3 2 6 5 3 5 1 2",
 
" 0 19  8 13 10 6 11 12 8 3 9 3 10 13 9 3 10 9 4 3 4 9 3 3 4 3 0 3 0 3 2 3 2 3 8 3 2 8 7 3 7 8 11 3 8 12 11 3 7 11 6 3 7 6 1 3 2 7 1 3 0 2 1 3 1 6 5 3 5 6 10 3 5 10 4 3 5 4 0 3 0 1 5",
 
" 0 17  8 12 11 8 4 5 1 6 9 3 11 12 9 3 11 9 10 3 10 9 6 3 10 6 7 3 7 6 2 3 2 6 1 3 2 1 0 3 0 1 5 3 0 5 4 3 0 4 3 3 2 0 3 3 7 2 3 3 3 4 8 3 7 3 8 3 10 7 8 3 8 11 10",
 
" 0 15  8 0 1 3 8 11 10 7 2 3 1 0 2 3 1 2 5 3 5 2 6 3 6 2 7 3 6 7 10 3 6 10 9 3 5 6 9 3 9 10 8 3 10 11 8 3 9 8 4 3 4 8 3 3 5 9 4 3 1 5 4 3 4 3 1",
 
" 0 15  8 8 4 1 0 2 5 9 11 3 8 11 10 3 10 11 9 3 10 9 5 3 10 5 6 3 8 10 6 3 6 5 2 3 6 2 7 3 8 6 7 3 4 8 7 3 7 2 3 3 3 2 0 3 4 7 3 3 1 4 3 3 3 0 1",
 
" 0 19  8 1 0 2 9 10 12 13 8 3 1 8 7 3 7 8 13 3 7 13 12 3 7 12 6 3 1 7 6 3 6 12 11 3 11 12 10 3 11 10 4 3 4 10 3 3 3 10 9 3 3 9 2 3 3 2 0 3 4 3 0 3 4 0 5 3 5 0 1 3 5 1 6 3 5 6 11 3 11 4 5", //  };

/* ADDED Oct 26, 1997:
};

 -- sort $Revision: 1.2 $ -- 
size = 17
recent size = 25
 -- arrange $Revision: 1.2 $ -- 
 -- arrangeInvariant $Revision: 1.2 $ -- 
//// combined size = 32
//// addition size = 15
public class graphUNKNOWN {
final static String data[] = {
*/

" 0 11  8 2 4 8 9 7 3 1 0 3 2 0 3 3 0 1 3 3 2 3 6 3 6 3 7 3 6 7 9 3 6 9 5 3 2 6 5 3 4 2 5 3 8 4 5 3 5 9 8",

" 0 13  8 8 6 1 2 3 7 9 10 3 8 10 7 3 10 9 7 3 8 7 4 3 4 7 3 3 4 3 0 3 0 3 2 3 0 2 1 3 0 1 5 3 4 0 5 3 8 4 5 3 6 8 5 3 5 1 6",

" 0 11  8 7 4 1 0 2 5 8 9 3 7 9 8 3 7 8 5 3 7 5 6 3 4 7 6 3 6 5 2 3 6 2 3 3 4 6 3 3 3 2 0 3 3 0 1 3 1 4 3",

" 0 11  8 0 2 6 7 8 9 5 1 3 0 1 5 3 0 5 4 3 4 5 8 3 5 9 8 3 4 8 7 3 4 7 3 3 0 4 3 3 2 0 3 3 6 2 3 3 3 7 6",

" 0 13  8 4 8 10 7 2 0 1 3 3 4 3 1 3 4 1 5 3 5 1 2 3 1 0 2 3 5 2 6 3 6 2 7 3 6 7 10 3 6 10 9 3 5 6 9 3 4 5 9 3 8 4 9 3 9 10 8",

" 0 11  8 6 9 8 3 4 0 1 5 3 6 5 1 3 6 1 2 3 2 1 0 3 2 0 3 3 0 4 3 3 2 3 7 3 7 3 8 3 6 2 7 3 9 6 7 3 7 8 9",

" 0 11  8 6 2 0 3 4 7 8 9 3 6 9 8 3 6 8 5 3 5 8 7 3 5 7 4 3 5 4 1 3 6 5 1 3 2 6 1 3 1 4 0 3 4 3 0 3 0 2 1",

" 0 11  8 4 0 1 2 5 9 7 8 3 4 8 7 3 4 7 6 3 6 7 9 3 6 9 5 3 6 5 3 3 4 6 3 3 0 4 3 3 3 5 2 3 0 3 2 3 2 1 0",

" 0 13  8 7 8 9 4 0 5 6 10 3 7 10 6 3 7 6 1 3 1 6 5 3 1 5 0 3 1 0 2 3 7 1 2 3 8 7 2 3 2 0 3 3 3 0 4 3 8 2 3 3 9 8 3 3 3 4 9",

" 0 13  8 6 10 8 3 4 0 5 9 3 6 9 5 3 6 5 1 3 1 5 0 3 1 0 2 3 2 0 3 3 0 4 3 3 2 3 8 3 2 8 7 3 1 2 7 3 6 1 7 3 10 6 7 3 7 8 10",

" 0 11  8 7 4 1 0 2 5 8 9 3 7 9 8 3 7 8 6 3 6 8 5 3 6 5 3 3 3 5 2 3 3 2 1 3 2 0 1 3 3 1 4 3 6 3 4 3 4 7 6",

" 0 13  8 4 8 9 10 7 2 0 3 3 4 3 0 3 4 0 1 3 1 0 2 3 1 2 6 3 6 2 7 3 6 7 9 3 7 10 9 3 6 9 5 3 5 9 8 3 1 6 5 3 4 1 5 3 5 8 4",

" 0 13  8 4 1 0 2 5 10 8 9 3 4 9 8 3 4 8 7 3 7 8 6 3 6 8 10 3 6 10 5 3 6 5 2 3 7 6 2 3 7 2 3 3 3 2 0 3 4 7 3 3 1 4 3 3 3 0 1",

" 0 13  8 3 6 10 9 5 2 0 1 3 3 1 4 3 4 1 2 3 1 0 2 3 4 2 5 3 4 5 8 3 3 4 8 3 8 5 9 3 8 9 7 3 3 8 7 3 6 3 7 3 10 6 7 3 7 9 10",

" 0 11  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 5 3 8 5 4 3 0 8 4 3 0 4 9 3 1 0 9 3 9 4 3 3 9 3 2 3 2 1 9",

/* Added Nov 4, 1997 :
};
 -- sort $Revision: 1.2 $ -- 
size = 32
recent size = 34
 -- arrange $Revision: 1.2 $ -- 
 -- arrangeInvariant $Revision: 1.2 $ -- 
//// combined size = 44
//// addition size = 12
public class graphUNKNOWN {
final static String data[] = {
*/

" 0 14  8 0 1 2 3 4 5 6 7 4 0 7 6 8 3 8 6 9 3 9 6 5 3 9 5 4 3 9 4 10 3 8 9 10 3 10 4 3 3 10 3 11 3 8 10 11 3 0 8 11 3 11 3 2 3 0 11 2 3 2 1 0",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 6 4 0 6 5 8 3 8 5 9 3 9 5 4 3 9 4 3 3 9 3 10 3 8 9 10 3 10 3 2 3 10 2 11 3 8 10 11 3 0 8 11 3 1 0 11 3 11 2 1",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 5 3 7 6 5 3 8 5 4 3 8 4 9 3 9 4 3 3 9 3 10 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 4 9 10 11 0 3 0 8 9",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 5 3 8 5 9 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 4 9 10 11 0 3 0 8 9",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 5 3 8 5 9 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 11 3 11 3 2 3 11 2 1 3 11 1 0 4 9 10 11 0 3 0 8 9",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 9 3 9 6 5 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 2 3 10 2 11 3 11 2 1 3 11 1 0 3 11 0 8 4 8 9 10 11",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 6 3 8 6 9 3 9 6 5 3 9 5 4 3 9 4 10 3 10 4 3 3 10 3 11 3 11 3 2 3 11 2 1 3 11 1 0 4 10 11 0 8 3 8 9 10",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 4 11 3 0 8 3 0 3 2 3 8 9 10 3 10 11 8 3 2 1 0",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 3 11 3 1 3 3 2 1 3 11 1 0 4 10 11 0 8 3 8 9 10",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 3 11 3 2 4 11 2 0 8 3 2 1 0 3 8 9 10 3 10 11 8",

" 0 16  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 5 3 9 5 10 3 10 5 4 3 10 4 11 3 11 4 3 3 11 3 2 3 11 2 12 3 12 2 1 3 12 1 0 3 12 0 8 4 10 11 12 8 3 8 9 10",

" 0 14  8 0 1 2 3 4 5 6 7 3 0 7 8 3 8 7 9 3 9 7 6 3 9 6 10 3 10 6 5 3 10 5 4 3 10 4 11 3 11 4 3 4 11 3 0 8 3 0 3 2 3 9 10 11 3 11 8 9 3 2 1 0"};

};