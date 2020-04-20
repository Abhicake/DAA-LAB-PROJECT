/*******************************************************************************
* Program: huffcode.c
* Purpose: To Huffman encode or decode a file. The program is incomplete because
* the student needs to code the functions in the file mytable.c.
* Name: Michael Heath, University of South Florida
* Date: 3/29/2000
* To Compile: gcc -o huffcode huffcode.c mytable.c imageio.c -lm
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "imageio.h"
#include "huffcode.h"

/*******************************************************************************
* This is a structure used to store the decoding tree.
*******************************************************************************/
struct DECODER{
   short int isleaf, value;
   struct DECODER *nextbit[2];
};

int VERBOSE=0;

static unsigned char queue[8];
static int numinqueue=0;

static unsigned char byte;
const unsigned char bmask=128;
static int availablebits=0;

void flush_queue();

int main(int argc, char *argv[])
{
   FILE *fp;
   unsigned char **image=NULL, **imageuc=NULL;
   int rows, cols, r, c, rowsuc, colsuc;
   unsigned long int original_length, compressed_length;
   char *task=NULL, *mode=NULL, *inputfilename=NULL, outputfilename[100];

   void get_commandline_parameters(int argc, char *argv[], char **task, char **mode,
       char **inputfilename, char *outputfilename);
   int write_huffman_image(char *filename, unsigned char **image, int rows, int cols);
   int read_huffman_image(char *filename, unsigned char ***image, int *rows, int *cols);
   int huffman_encode_datafile(char *inputfilename, char *outputfilename);
   int huffman_decode_datafile(char *inputfilename, char *outputfilename);

   /****************************************************************************
   * Check the list of command line arguments and offer help to the user if
   * they only enter the name of the program with no other arguments.
   ****************************************************************************/
   get_commandline_parameters(argc, argv, &task, &mode, &inputfilename, outputfilename);

   if(VERBOSE){
      printf("inputfilename = %s\n", inputfilename);
      printf("task = %s\n", task);
      printf("mode = %s\n", mode);
      printf("outputfilename = %s\n", outputfilename);
   }

   /****************************************************************************
   * Compress the data if that is what we are supposed to do.
   ****************************************************************************/
   if(strcmp(task, "COMPRESS") == 0){

      /*************************************************************************
      * If we are to transform the file and compress it, do that.
      *************************************************************************/
      if(strcmp(mode, "TRANSFORM") == 0){
         if(read_pgm_image(inputfilename, &image, &rows, &cols) == 0) exit(1);
         write_huffman_image(outputfilename, image, rows, cols);
	 free_image(image, rows);
      }

      /*************************************************************************
      * If we are to Huffman encode the raw file data, do it.
      *************************************************************************/
      if(strcmp(mode, "DATA") == 0){
         huffman_encode_datafile(inputfilename, outputfilename);
      }
   }

   /****************************************************************************
   * Uncompress the data if that is what we are supposed to do.
   ****************************************************************************/
   if(strcmp(task, "UNCOMPRESS") == 0){

      /*************************************************************************
      * If the compressed file was transformed, decompress the image.
      *************************************************************************/
      if(strcmp(mode, "TRANSFORM") == 0){
         read_huffman_image(inputfilename, &image, &rows, &cols);
         write_pgm_image(outputfilename, image, rows, cols, (char *)NULL, 255);
      }

      /*************************************************************************
      * If the compressed file was not transformed, decompress the data.
      *************************************************************************/
      if(strcmp(mode, "DATA") == 0){
         huffman_decode_datafile(inputfilename, outputfilename);
      }
   }
}

/*******************************************************************************
* Function: write_huffman_image
* Purpose: This function will save an image in a file. The image is compressed
* by Huffman coding. A "transformation" is used to attempt to get a better
* compression by reducing the entropy (disorder) in the image. Nearby pixel
* values are used to predict the value of the current pixel, and the difference
* between the predicted value and the actual value is what gets encoded.
*
* A map of the image showing which         This is a local window that shows the
* mathematical equation is used to         variables to use in the prediction
* predict the pixel in the image.          equation. Place the window with "a"
*                                          laying on the pixel for which the
*      0 1 2 3 4 5 6 7 8                   estimate is to be computed.
*     +-+-+-+-+-+-+-+-+-+
*   0 |A|B|B|B|B|B|B|B|B|                       +-+-+-+        
*     +-+-+-+-+-+-+-+-+-+                       |e|d|c|        
*   1 |C|D|D|D|D|D|D|D|E|                       +-+-+
*     +-+-+-+-+-+-+-+-+-+                       |b|a|
*   2 |C|D|D|D|D|D|D|D|E|                       +-+-+
*     +-+-+-+-+-+-+-+-+-+
*   3 |C|D|D|D|D|D|D|D|E|
*     +-+-+-+-+-+-+-+-+-+
*
* At each pixel position in the image we calculate an estimate (a_est) of the
* pixels value using other pixels that would have already been seen in a raster
* scan of the image up to the position of the current pixel in the image. A
* difference between the current pixel and the estimate is what will be encoded
* with the Huffman code. This "transformed image" will often be more
* "compressable" than the original image was (when Huffman encoding is used).
*
*  A: a_est = a
*  B: a_est = b
*  C: a_est = d
*  D: a_est = floor((b+c+d+e)/4)
*  E: a_est = floor((b+d+e)/3)
*
* Name: Michael Heath, University of South Florida 
*******************************************************************************/
int write_huffman_image(char *filename, unsigned char **image, int rows, int cols)
{
   FILE *fp;
   struct HUFFTABLE *table=NULL;
   int i, k, r, c, maxval, minval, numval;
   unsigned long int headerbytes, tablebytes;
   short int **imdiff=NULL;

   void form_table(struct HUFFTABLE *table, int num);
   void write_huffman_table(struct HUFFTABLE *table, int min, int max, FILE *fp);
   void write_bits(int numbits, unsigned char *bitstring, FILE *fp);

   /****************************************************************************
   * "Transform" the image. The input image is unsigned char with values in
   * the range 0->255. The "transformed" image could have values in the
   * range -255->255. Therefore we will use an unsigned short integer image
   * to store the "transformed" image.
   ****************************************************************************/
   if((imdiff = (short int **) calloc(rows, sizeof(short int *))) == NULL){
      fprintf(stderr, "Calloc error in write_huffman_image().\n");
      return(0);
   }
   for(r=0;r<rows;r++){
      if((imdiff[r] = (short int *) calloc(cols, sizeof(short int))) == NULL){
         fprintf(stderr, "Calloc error in write_huffman_image().\n");
         return(0);
      }
   }

   /****************************************************************************
   * Compute the "transformed" image.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      if(r == 0){
         imdiff[r][0] = image[r][0];
         for(c=1;c<cols;c++){
            imdiff[r][c] = (short int)image[r][c] - (short int)image[r][c-1];
         }
      }
      else{
         c = 0;
         imdiff[r][c] = (short int)image[r][c] - (short int)image[r-1][c];
         for(c=1;c<(cols-1);c++)
            imdiff[r][c] = (short int)image[r][c] -
               (short int)(((int)image[r][c-1] + (int)image[r-1][c-1] +
                            (int)image[r-1][c] + (int)image[r-1][c+1]) / 4);
         c = cols-1;
         imdiff[r][c] = (short int)image[r][c] -
            (short int)(((int)image[r][c-1] + (int)image[r-1][c-1] +
                         (int)image[r-1][c]) / 3);
      }
   }

   /****************************************************************************
   * Determine the minimum and maximum value in the imdiff image and allocate
   * a Huffman code table large enough to store all of the entries.
   ****************************************************************************/
   maxval = minval = imdiff[0][0];
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         if(imdiff[r][c] > maxval) maxval = imdiff[r][c];
         else if(imdiff[r][c] < minval) minval = imdiff[r][c];
      }
   }
   if((table = (struct HUFFTABLE *) calloc(maxval-minval+1, sizeof(struct HUFFTABLE))) == NULL){
      fprintf(stderr, "Calloc error in write_huffman_image()!\n");
      return(0);
   }

   /****************************************************************************
   * Open a file in which we will save the Huffman encoded image.
   ****************************************************************************/
   if((fp = fopen(filename, "wb")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing!\n", filename);
      free(table);
      return(0);
   }

   /****************************************************************************
   * Compute the frequencies of pixels and store them in a structure with
   * other fields of data. This will end up being the Huffman code table.
   ****************************************************************************/
   for(i=0;i<=(maxval-minval);i++){
      table[i].value = minval+i;
      /*
      table[i].maxbits = 8;
      table[i].bitstring = (unsigned char *) calloc(8, sizeof(unsigned char));
      */
      table[i].maxbits = 0;
      table[i].bitstring = (unsigned char *)NULL;
   }
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         table[imdiff[r][c]-minval].count += 1;
      }
   }

   /****************************************************************************
   * Forn the Huffman code table.
   ****************************************************************************/
   form_table(table, maxval-minval);

   /****************************************************************************
   * Just for fun, compute the average bit rate for the encoded image (less the
   * header and Huffman table) and the entropy. 
   ****************************************************************************/
   if(VERBOSE){
      double prob, entropy;
      unsigned long int sum, totalcount;

      sum = 0;
      totalcount = 0;
      for(i=0;i<=(maxval-minval);i++){
         sum += table[i].count * table[i].bits;
         totalcount += table[i].count;
      }
      entropy = 0.0;
      for(i=0;i<=(maxval-minval);i++){
         if(table[i].count != 0){
            prob = (double)table[i].count / (double)(totalcount);
            entropy -= prob * log10(prob)/log10(2.0);
         }
      }
      printf("The entropy of the image data is %f.\n", entropy);
      printf("The avgerage bit rate of the image data is %f\n", (float)sum / (float)totalcount);
   }

   /****************************************************************************
   * Write the header information to the file.
   ****************************************************************************/
   for(i=0,numval=0;i<=(maxval-minval);i++) if(table[i].count != 0) numval++;
   fprintf(fp, "HFT\n%d %d\n%d %d\n%d\n", cols, rows, minval, maxval, numval);
   headerbytes = ftell(fp);
   write_huffman_table(table, minval, maxval, fp);
   tablebytes = ftell(fp) - headerbytes;

   if(VERBOSE){
      printf("The header takes %lu bytes and the Huffman table takes %lu bytes.\n",
         headerbytes, tablebytes);
   }

   /****************************************************************************
   * Encode the image, writing it to the file as we go.
   ****************************************************************************/
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         write_bits(table[imdiff[r][c]-minval].bits, table[imdiff[r][c]-minval].bitstring, fp);
      }
   }
   if(numinqueue != 0) flush_queue(fp);

   /****************************************************************************
   * Free up the memory for the difference image.
   ****************************************************************************/
   for(r=0;r<rows;r++) free(imdiff[r]);
   free(imdiff);

   free(table);
   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: write_huffman_table
* Purpose: To write the Huffman table out to a file.
* Name: Michael Heath, University of South Florida
* Date: 1/11/2000
*******************************************************************************/
void write_huffman_table(struct HUFFTABLE *table, int min, int max, FILE *fp)
{
   int i;
   unsigned char achar;

   void write_bits(int numbits, unsigned char *bitstring, FILE *fp);
   void flush_queue(FILE *fp);

   for(i=0;i<=(max-min);i++){
      if(table[i].count != 0){
         /* fwrite(&(table[i].value), sizeof(short int), 1, fp); */
         achar = table[i].value >> 8;
         fputc(achar, fp);
         achar = table[i].value % 256;
         fputc(achar, fp);
         fputc((unsigned char)table[i].bits, fp);
         write_bits(table[i].bits, table[i].bitstring, fp);
         flush_queue(fp);
      }
   }
}

/*******************************************************************************
* Function: write_bits
* Purpose: To write the elements in an array of unsigned characters out to
* a file as a bit stream.
* Name: Michael Heath, University of South Florida
* Date: 1/11/2000
*******************************************************************************/
void write_bits(int numbits, unsigned char *bitstring, FILE *fp)
{
   int i;

   void flush_queue(FILE *fp);

   for(i=0;i<numbits;i++){
      queue[numinqueue] = bitstring[i];
      numinqueue++;
      if(numinqueue == 8) flush_queue(fp);
   }
}

/*******************************************************************************
* Function: flush_queue
* Purpose: To flush the values in the queue array out to a file.
* Name: Michael Heath, University of South Florida
* Date: 1/11/2000
*******************************************************************************/
void flush_queue(FILE *fp)
{
   unsigned char achar=0;

   if(numinqueue == 0) return;

   numinqueue = 0;

   if(queue[0] == 1) achar += 128;
   if(queue[1] == 1) achar += 64;
   if(queue[2] == 1) achar += 32;
   if(queue[3] == 1) achar += 16;
   if(queue[4] == 1) achar += 8;
   if(queue[5] == 1) achar += 4;
   if(queue[6] == 1) achar += 2;
   if(queue[7] == 1) achar += 1;

   fputc(achar, fp);

   memset(queue, 0, 8);
}

/*******************************************************************************
* Function: get_huffman_pixel
* Purpose: The function will return the next decoded pixel value in the
* image. It is a recursive function.
* Name: Michael Heath, University of South Florida
* Date: 1/13/2000
*******************************************************************************/
int get_huffman_pixel(FILE *fp, struct DECODER *decoder)
{
   unsigned char bval;

   if(decoder->isleaf == 1) return(decoder->value);

   if(availablebits == 0){
      byte = fgetc(fp);
      availablebits = 8;
   }

   if(byte & bmask) bval = 1;
   else bval = 0;

   byte = byte<<1;
   availablebits--;

   if(bval == 0) return(get_huffman_pixel(fp, decoder->nextbit[0]));
   else return(get_huffman_pixel(fp, decoder->nextbit[1]));
}

/*******************************************************************************
* Function: make_decode_tree
* Purpose: We have the Huffman code table that provides the bit code for each
* pixel value, but we don't have an efficient way to find the pixel value
* for a bit code. This function forms a tree which can be traversed with the
* bit from a bit stream to produce the decoded pixel values. Most all of the
* work of doing this is accomplished by the recursive function addleaf().
* Name: Michael Heath, University of South Florida
* Date: 1/12/2000
*******************************************************************************/
void make_decode_tree(struct HUFFTABLE *table, int tablelen, struct DECODER **decodertree)
{
   int i;
   struct DECODER *decoder=NULL;

   void addleaf(struct DECODER *decoder, int len, unsigned char *arr, int value);

   decoder = (struct DECODER *) calloc(1, sizeof(struct DECODER));
   *decodertree = decoder;

   for(i=0;i<tablelen;i++){
      addleaf(decoder, table[i].bits, table[i].bitstring, table[i].value);
   }
}

/*******************************************************************************
* Function: addleaf
* Purpose: This function adds a leaf node to a binary tree. The leaf node
* represents a pixel value. The path from the root of the tree to the leaf
* is the Huffman bit code for the pixel value.
* Name: Michael Heath, University of South Florida
* Date: 1/12/2000
*******************************************************************************/
void addleaf(struct DECODER *decoder, int len, unsigned char *arr, int value)
{
   if(len == 0){
      decoder->isleaf = 1;
      decoder->value = value;
      return;
   }

   if(decoder->nextbit[*arr] == NULL) decoder->nextbit[*arr] = calloc(1, sizeof(struct DECODER));
   addleaf(decoder->nextbit[*arr], len-1, arr+1, value);
}

/*******************************************************************************
* Function: show_decoder_tree
* Purpose: This function can display the contents of the decoder tree. It
* performs a breadth first search on the tree to print out the paired
* Huffman bit codes and pixel values. This is a recursive function. The
* first call to the function must initialize k to zero and provide an
* array showarr of length maxval+1.
* Name: Michael Heath, University of South Florida
* Date: 1/13/2000
*******************************************************************************/
void show_decoder_tree(struct DECODER *decoder, int *k, char *showarr)
{
   int i;

   if(decoder == NULL) return;

   if(decoder->nextbit[0] != NULL){
      showarr[*k] = 0;
      (*k)++;
      show_decoder_tree(decoder->nextbit[0], k, showarr);
      (*k)--;
   }

   if(decoder->nextbit[1] != NULL){
      showarr[*k] = 1;
      (*k)++;
      show_decoder_tree(decoder->nextbit[1], k, showarr);
      (*k)--;
   }

   if(decoder->isleaf == 1){
      for(i=0;i<(*k);i++) printf("%d", showarr[i]);
      printf("<%d>", decoder->value);
      printf("\n");
      return;
   }
}

/*******************************************************************************
* Function: read_huffman_image
* Purpose: This function will read in an image that is was stored to a file
* using the function write_huffman_image().
* Name: Michael Heath, University of South Florida
* Date: 1/12/2000
*******************************************************************************/
int read_huffman_image(char *filename, unsigned char ***image, int *rows, int *cols)
{
   FILE *fp;
   struct HUFFTABLE *table = NULL;
   int i, k, r, c, minval, maxval, numval;
   unsigned char **im=NULL;
   struct DECODER *decoder=NULL;
   char *showarr=NULL;
   int delta;

   void read_huffman_table(struct HUFFTABLE *table, int numval, FILE *fp);
   void make_decode_tree(struct HUFFTABLE *table, int tablelen, struct DECODER **decodertree);
   void show_decoder_tree(struct DECODER *decoder, int *k, char *showarr);
   int get_huffman_pixel(FILE *fp, struct DECODER *decoder);

   /****************************************************************************
   * Open the file containing the Huffman encoded image.
   ****************************************************************************/
   if((fp = fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s for reading!\n", filename);
      return(0);
   }

   /****************************************************************************
   * Read in the image header and the huffman code table.
   ****************************************************************************/
   fscanf(fp, "HFT\n%d %d\n%d %d\n%d\n", cols, rows, &minval, &maxval, &numval);

   /****************************************************************************
   * Allocate an array to store the Huffman code table we want to read from the
   * file. Then read it in from the file.
   ****************************************************************************/
   if((table = (struct HUFFTABLE *) calloc(numval, sizeof(struct HUFFTABLE))) == NULL){
      fprintf(stderr, "Calloc error in read_huffman_image()!\n");
      return(0);
   }
   read_huffman_table(table, numval, fp);
   make_decode_tree(table, numval, &decoder);

   /****************************************************************************
   * Allocate the image.
   ****************************************************************************/
   if((im = allocate_image(*rows, *cols)) == NULL) exit(1);
   *image = im;

   /****************************************************************************
   * Read in the bit stream, decomressing it is read from the file. Un-transform
   * the image as we read it.
   ****************************************************************************/
   availablebits = 0;

   for(r=0;r<(*rows);r++){
      if(r == 0){
         im[r][0] = (unsigned char)get_huffman_pixel(fp, decoder);
         for(c=1;c<(*cols);c++){
            delta = get_huffman_pixel(fp, decoder);
            im[r][c] = (unsigned char)((int)im[r][c-1] + delta);
         }
      }
      else{
         c = 0;
         im[r][c] = (unsigned char)((int)im[r-1][c] + get_huffman_pixel(fp, decoder));
         for(c=1;c<((*cols)-1);c++){
            im[r][c] = (unsigned char)((((int)im[r][c-1] + (int)im[r-1][c-1] + (int)im[r-1][c] + (int)im[r-1][c+1]) / 4)
			    + get_huffman_pixel(fp, decoder));
         }
         c = (*cols) -1;
         im[r][c] = (unsigned char)((((int)im[r][c-1] + (int)im[r-1][c-1] + (int)im[r-1][c]) / 3)
			    + get_huffman_pixel(fp, decoder));
      }
   }

   free(table);
   fclose(fp);
   return(1);
}

/*******************************************************************************
* Function: read_huffman_table 
* Purpose: To read the huffman table from a file.
* Name: Michael Heath
* Date: 1/12/2000
*******************************************************************************/
void read_huffman_table(struct HUFFTABLE *table, int numval, FILE *fp)
{
   int i, b, numbits;
   short int value;
   unsigned char bval;
   unsigned short int val;

   availablebits = 0;
   for(i=0;i<numval;i++){
      /* fread(&value, sizeof(short int), 1, fp); */
      value = fgetc(fp);
      value = value<<8;
      value += fgetc(fp);
      numbits = fgetc(fp);
      table[i].value = value;
      table[i].bits = numbits;

      if(numbits != 0){

         table[i].bitstring = (unsigned char *) calloc(numbits, 1);

         b = 0;
         availablebits = 0;
         while(b < numbits){
            if(availablebits == 0){
               byte = fgetc(fp);
               availablebits = 8;
            }

            if(byte & bmask) bval = 1;
            else bval = 0;
            byte = byte<<1;
            availablebits--;

            table[i].bitstring[b] = bval;
            b++;
         }
      }
   }
}

/*******************************************************************************
* Function: print_help
* Purpose: This function prints out help on using this program.
* Name: Michael Heath, University of South Florida
* Date: 1/14/2000
*******************************************************************************/
void print_help(char *string)
{
   fprintf(stderr, "\n\n********************************************************************************\n");
   fprintf(stderr, "This program implements a lossless encoding using Huffman compression. The\n");
   fprintf(stderr, "program can compress or decompress a file. The compression can be applied to\n");
   fprintf(stderr, "a file, interpreting it as raw data by using the -data option. If the input\n");
   fprintf(stderr, "file is an image stored in the RAW PGM format it can be compressed by\n");
   fprintf(stderr, "transforming the file to reduce the redundancy before Huffman encoding it\n");
   fprintf(stderr, "by using the option -pgm. The default mode is to examine the contents of\n");
   fprintf(stderr, "the input file. If the file appears to contain compressed data, it will\n");
   fprintf(stderr, "be uncompressed. If the file does not appear to contain compressed data, it\n");
   fprintf(stderr, "will be compressed in the -data mode. If the data is being compressed, it\n");
   fprintf(stderr, "will be written to a file whose name is the input filename followed by an\n");
   fprintf(stderr, "appended .hf extension. If the file is being decompressed, the\n");
   fprintf(stderr, "decompressed data will be written to a file with the input filename with\n");
   fprintf(stderr, "an extension of .hfu appended to it.\n");
   fprintf(stderr, "********************************************************************************\n\n");
   fprintf(stderr, "<USAGE> %s [-pgm | -data | -uc] filename\n\n", string);
}

/*******************************************************************************
* Function: get_commandline_parameters
* Purpose: This function gets the command line parameters and checks to make
* sure they are valid. It also helps determine what mode to run the code in
* if it is not explicitly specified by the user.
* Name: Michael Heath, University of South Florida
* Date: 1/14/2000
*******************************************************************************/
void get_commandline_parameters(int argc, char *argv[], char **task, char **mode,
    char **inputfilename, char *outputfilename)
{
   FILE *fp;
   int i;
   char char1='\0', char2='\0', char3='\0';

   void print_help(char *string);

   if(argc == 1){
      print_help(argv[0]);
      exit(1);
   }
   for(i=1;i<argc;i++){
      if(strcmp(argv[i], "-v") == 0){
         VERBOSE = 1;
      }
      else if(argv[i][0] == '-'){
         if(!((*task == NULL) || (*mode == NULL))){
            fprintf(stderr, "Error! You can not specify more than one option!\n");
            print_help(argv[0]);
            exit(1);
         }
         if(strcmp(argv[i], "-pgm") == 0){ *task = "COMPRESS"; *mode = "TRANSFORM"; }
         if(strcmp(argv[i], "-data") == 0){ *task = "COMPRESS"; *mode = "DATA"; }
         if(strcmp(argv[i], "-uc") == 0) *task = "UNCOMPRESS";
      }
      else{
         if(*inputfilename != NULL){
            fprintf(stderr, "Error! You can not specify more than filename!\n");
            print_help(argv[0]);
            exit(1);
         }
         *inputfilename = argv[i];
      }
   }

   /****************************************************************************
   * Check to see that the user input a filename.
   ****************************************************************************/
   if(*inputfilename == NULL){
      fprintf(stderr, "Error! You did not specify an input filename!\n");
      print_help(argv[0]);
      exit(1);
   }

   /****************************************************************************
   * If the user did not specify a mode or they specified, check the file to see
   * if if looks like compressed data or compressed transformed data.
   ****************************************************************************/
   if((*task == NULL) || (strcmp(*task, "UNCOMPRESS") == 0)){
      if((fp = fopen(*inputfilename, "rb")) == NULL){
         fprintf(stderr, "Error opening the file %s for reading!\n\n", *inputfilename);
         exit(1);
      }
      char1 = fgetc(fp);
      char2 = fgetc(fp);
      char3 = fgetc(fp);

      if((char1 == 'H') && (char2 == 'F')){
         if(char3 == 'D'){
            *task = "UNCOMPRESS";
            *mode = "DATA";
         }
         else if(char3 == 'T'){
            *task = "UNCOMPRESS";
            *mode = "TRANSFORM";
         }
         else{
            fprintf(stderr, "You asked to decomress the file and it does not appear to be compressed!\n");
            exit(1);
         }
      }
      else{
         if(*task == NULL){
            *task = "COMPRESS";
            *mode = "DATA";
         }
         else{
            fprintf(stderr, "You asked to decomress the file and it does not appear to be compressed!\n");
            exit(1);
         }
      }
      fclose(fp);
   }
   else{
      *task = "COMPRESS";
      if(*mode == NULL) *mode = "DATA";
   }

   if(strcmp(*task, "COMPRESS") == 0) sprintf(outputfilename, "%s.hf", *inputfilename);
   if(strcmp(*task, "UNCOMPRESS") == 0) sprintf(outputfilename, "%s.hfu", *inputfilename);
}

/*******************************************************************************
* Function: huffman_encode_datafile
* Purpose: To Huffman encode the data in a file.
* Name: Michael Heath, University of South Florida
* Date: 1/14/2000
*******************************************************************************/
int huffman_encode_datafile(char *inputfilename, char *outputfilename)
{
   FILE *fpin, *fpout;
   struct HUFFTABLE *table=NULL;
   int i, k, r, c, maxval, minval=0, numval;
   unsigned long int headerbytes, tablebytes, totalbytes=0;
   short int byteval;

   void form_table(struct HUFFTABLE *table, int num);
   void write_huffman_table(struct HUFFTABLE *table, int min, int max, FILE *fp);
   void write_bits(int numbits, unsigned char *bitstring, FILE *fp);


   /****************************************************************************
   * Allocate a Huffman code table large enough to store all of the entries.
   ****************************************************************************/
   maxval = 255; /* We know this because we are interpreting the file as bytes. */
   if((table = (struct HUFFTABLE *) calloc(maxval+1, sizeof(struct HUFFTABLE))) == NULL){
      fprintf(stderr, "Calloc error in huffman_encode_datafile()!\n");
      return(0);
   }

   /****************************************************************************
   * Compute the frequencies of pixels in the input file.
   ****************************************************************************/
   if((fpin = fopen(inputfilename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s for reading!\n", inputfilename);
      free(table);
      return(0);
   }
   byteval = fgetc(fpin);
   while(byteval != EOF){
      totalbytes++;
      table[byteval].count++;
      byteval = fgetc(fpin);
   }
   fclose(fpin);
   for(i=0;i<=maxval;i++){
      table[i].value = minval+i;
      table[i].maxbits = 8;
      table[i].bitstring = (unsigned char *) calloc(8, sizeof(unsigned char));
   }

   /****************************************************************************
   * Open a file in which we will save the Huffman encoded image.
   ****************************************************************************/
   if((fpout = fopen(outputfilename, "wb")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing!\n", outputfilename);
      free(table);
      return(0);
   }

   /****************************************************************************
   * Form the Huffman code table.
   ****************************************************************************/
   form_table(table, maxval);

   /****************************************************************************
   * Just for fun, compute the average bit rate for the encoded data (less the
   * header and Huffman table) and the entropy. 
   ****************************************************************************/
   if(VERBOSE){
      double prob, entropy;
      unsigned long int sum, totalcount;

      sum = 0;
      totalcount = 0;
      for(i=0;i<=maxval;i++){
         sum += table[i].count * table[i].bits;
         totalcount += table[i].count;
      }
      entropy = 0.0;
      for(i=0;i<=maxval;i++){
         if(table[i].count != 0){
            prob = (double)table[i].count / (double)(totalcount);
            entropy -= prob * log10(prob)/log10(2.0);
         }
      }
      printf("The entropy of the data is %f.\n", entropy);
      printf("The avgerage bit rate of the data is %f\n", (float)sum / (float)totalcount);
   }

   /****************************************************************************
   * Write the header information to the file.
   ****************************************************************************/
   for(i=0,numval=0;i<=maxval;i++) if(table[i].count != 0) numval++;
   fprintf(fpout, "HFD\n%lu\n%d %d\n%d\n", totalbytes, minval, maxval, numval);
   headerbytes = ftell(fpout);
   write_huffman_table(table, minval, maxval, fpout);
   tablebytes = ftell(fpout) - headerbytes;

   if(VERBOSE){
      printf("The header takes %lu bytes and the Huffman table takes %lu bytes.\n",
         headerbytes, tablebytes);
   }

   /****************************************************************************
   * Read the data, encode it and write it to a file.
   ****************************************************************************/
   if((fpin = fopen(inputfilename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s for reading!\n", inputfilename);
      free(table);
      return(0);
   }
   byteval = fgetc(fpin);
   while(byteval != EOF){
      write_bits(table[byteval].bits, table[byteval].bitstring, fpout);
      byteval = fgetc(fpin);
   }
   if(numinqueue != 0) flush_queue(fpout);

   fclose(fpin);
   fclose(fpout);

   free(table);
   return(1);
}

/*******************************************************************************
* Function: huffman_decode_datafile
* Purpose: To decode a file that was Huffman encoded without transformation.
* Name: Michael Heath, University of South Florida
* Date: 1/14/2000
*******************************************************************************/
int huffman_decode_datafile(char *inputfilename, char *outputfilename)
{
   FILE *fpin, *fpout;
   struct HUFFTABLE *table = NULL;
   int i, k, r, c, minval, maxval, numval;
   struct DECODER *decoder=NULL;
   char *showarr=NULL;
   unsigned long int totalbytes=0;
   unsigned char abyte;

   void read_huffman_table(struct HUFFTABLE *table, int numval, FILE *fp);
   void make_decode_tree(struct HUFFTABLE *table, int tablelen, struct DECODER **decodertree);
   void show_decoder_tree(struct DECODER *decoder, int *k, char *showarr);
   int get_huffman_pixel(FILE *fp, struct DECODER *decoder);

   /****************************************************************************
   * Open the file containing the Huffman encoded data.
   ****************************************************************************/
   if((fpin = fopen(inputfilename, "rb")) == NULL){
      fprintf(stderr, "Error opening the file %s for reading!\n", inputfilename);
      return(0);
   }

   /****************************************************************************
   * Read in the image header and the huffman code table.
   ****************************************************************************/
   fscanf(fpin, "HFD\n%lu\n%d %d\n%d\n", &totalbytes, &minval, &maxval, &numval);

   /****************************************************************************
   * Allocate an array to store the Huffman code table we want to read from the
   * file. Then read it in from the file.
   ****************************************************************************/
   if((table = (struct HUFFTABLE *) calloc(numval, sizeof(struct HUFFTABLE))) == NULL){
      fprintf(stderr, "Calloc error in read_huffman_image()!\n");
      return(0);
   }
   read_huffman_table(table, numval, fpin);
   make_decode_tree(table, numval, &decoder);

   /****************************************************************************
   * Open the file to write the data to.
   ****************************************************************************/
   if((fpout = fopen(outputfilename, "wb")) == NULL){
      fprintf(stderr, "Error opening the file %s for writing!\n", outputfilename);
      exit(1);
   }

   /****************************************************************************
   * Read in the bit stream, decomressing it is read from the file. 
   ****************************************************************************/
   availablebits = 0;
   for(i=0;i<totalbytes;i++){
      abyte = (unsigned char)get_huffman_pixel(fpin, decoder);
      fputc(abyte, fpout);
   }

   free(table);
   fclose(fpin);
   fclose(fpout);
   return(1);
}
