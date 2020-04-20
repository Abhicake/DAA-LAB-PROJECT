/*******************************************************************************
* File: huffcode.h
* Purpose: To share the HUFFTABLE data structure with the code that students
* must write to create the Huffman code table.
* Name: Michael Heath, University of South Florida
* Date: 3/29/2000
*******************************************************************************/

/*******************************************************************************
* These is a structure used to store the Huffman code table.
*
* value     = The pixel value (or "transformed" pixel value) it is set in huffcode.c
*
* count     = The number of occurances of the associated "value" and it is set
*             in huffcode.c.
*
* bits      = The number of bits for the Huffman code for the associated
*             value. It must be set in the function form_table in mytable.c.
*
* maxbits   = A variable you can use to keep track of the length of the array
*             bitstring after you allocate it. You need not use this variable
*             if you do not want to. It is in the structure, but my code
*             will not care what value you place in this variable.
*
* bitstring = This is an unsigned character array that is used to store the
*             bit string of the Huffman code for the associated value. It must
*             be set in the function form_table in mytable.c.
*             An unsigned character array is used instead of actually setting
*             bits. Use (unsigned char)0 to represent the bit value of zero
*             and (unsigned char)1 to represent the bit value of one. Thus a
*             a bit code of 011 would be stored in an unsigned char array
*             with a length of at least 3 (THAT IS ALLOCATED (by the student)
*             IN form_table), bits would be set to 3 and the first
*             three elements in the array would have the values
*             (unsigned char)0, (unsigned char)1 and (unsigned char)1
*             respectively.
*******************************************************************************/
struct HUFFTABLE{
   short int value;
   long int count;
   int bits;
   int maxbits;
   unsigned char *bitstring;
};
