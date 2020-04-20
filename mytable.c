/*******************************************************************************
* File: mytable.c
* Purpose: To supply the function form_table() to the program.
* Name: 
* Date: 
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "huffcode.h"

extern int VERBOSE;

/*******************************************************************************
* Function: form_table
* Purpose: To calculate the Huffman code for pixel frequencies stored in the
* array of structures named table. The array has a length of num. The value
* and count data fields are already assigned values when the function is called.
* You should write the body of this function. To do so they need to
* implement an algorithm to find a Huffman code table for given frequencies,
* which are stored as the counts fields in the array named table). You
* should set the bits field to be the number of bits of the Huffman code, and
* fill the elements of the unsigned character array "bitstring" with the
* values (unsigned char)0 and (unsigned char)1 to store the Huffman code.
* Please note that the bitstring field is an unsigned char pointer and you
* should allocate memory for the array bitstring for each element of the
* table (use calloc or malloc).
*******************************************************************************/


#define MAX_TREE_HT 9999

// A Huffman tree node 
struct MinHeapNode { 

	// One of the input intacters 
	short int data; 

	// Frequency of the intacter 
	long int freq; 

	struct HUFFTABLE *add ;
	// Left and right child of this node 
	struct MinHeapNode *left, *right; 
}; 

// A Min Heap: Collection of 
// min-heap (or Huffman tree) nodes 
struct MinHeap { 

	// Current size of min heap 
	unsigned size; 

	// capacity of min heap 
	unsigned capacity; 

	// Array of minheap node pointers 
	struct MinHeapNode** array; 
}; 

// A utility function allocate a new 
// min heap node with given intacter 
// and frequency of the intacter 
struct MinHeapNode* newNode(short int data, long int freq,struct HUFFTABLE *address) 
{ 
	struct MinHeapNode* temp = (struct MinHeapNode*)malloc(sizeof(struct MinHeapNode)); 

	temp->left = temp->right = NULL; 
	temp->data = data; 
	temp->freq = freq; 
	temp->add = address;

	return temp; 
} 

// A utility function to create 
// a min heap of given capacity 
struct MinHeap* createMinHeap(unsigned capacity) 

{ 

	struct MinHeap* minHeap = (struct MinHeap*)malloc(sizeof(struct MinHeap)); 

	// current size is 0 
	minHeap->size = 0; 

	minHeap->capacity = capacity; 

	minHeap->array = (struct MinHeapNode**)malloc(minHeap->capacity * sizeof(struct MinHeapNode*)); 

	return minHeap; 
} 

// A utility function to 
// swap two min heap nodes 
void swapMinHeapNode(struct MinHeapNode** a, 
					struct MinHeapNode** b) 

{ 

	struct MinHeapNode* t = *a; 
	*a = *b; 
	*b = t; 
} 

// The standard minHeapify function. 
void minHeapify(struct MinHeap* minHeap, int idx) 

{ 

	int smallest = idx; 

	int left = 2 * idx + 1; 
	int right = 2 * idx + 2; 

	if (left < minHeap->size && minHeap->array[left]-> freq < minHeap->array[smallest]->freq) 
		smallest = left; 

	if (right < minHeap->size && minHeap->array[right]-> freq < minHeap->array[smallest]->freq) 
		smallest = right; 

	if (smallest != idx) { 
		swapMinHeapNode(&minHeap->array[smallest], 
						&minHeap->array[idx]); 
		minHeapify(minHeap, smallest); 
	} 
} 

// A utility function to check 
// if size of heap is 1 or not 
int isSizeOne(struct MinHeap* minHeap) 
{ 

	return (minHeap->size == 1); 
} 

// A standard function to extract 
// minimum value node from heap 
struct MinHeapNode* extractMin(struct MinHeap* minHeap) 

{ 

	struct MinHeapNode* temp = minHeap->array[0]; 
	minHeap->array[0] 
		= minHeap->array[minHeap->size - 1]; 

	--minHeap->size; 
	minHeapify(minHeap, 0); 
	
	return temp; 
} 

// A utility function to insert 
// a new node to Min Heap 
void insertMinHeap(struct MinHeap* minHeap, 
				struct MinHeapNode* minHeapNode) 

{ 

	++minHeap->size; 
	int i = minHeap->size - 1; 

	while (i && minHeapNode->freq < minHeap->array[(i - 1) / 2]->freq) { 

		minHeap->array[i] = minHeap->array[(i - 1) / 2]; 
		i = (i - 1) / 2; 
	} 

	minHeap->array[i] = minHeapNode; 
} 

// A standard function to build min heap 
void buildMinHeap(struct MinHeap* minHeap) 
{ 
		// for (int i = 0; i < 6; ++i) 
		// printf("3h: %d %ld %p\n", minHeap->array[i]->data, minHeap->array[i]->freq, minHeap->array[i]->add);
	
	int n = minHeap->size - 1; 
	int i; 
	for (i = (n - 1) / 2; i >= 0; --i) 
		minHeapify(minHeap, i); 

	// for (int i = 0; i < 6; ++i) 
	// 	printf("4h: %d %ld %p\n", minHeap->array[i]->data, minHeap->array[i]->freq, minHeap->array[i]->add);
	
} 


// Utility function to check if this node is leaf 
int isLeaf(struct MinHeapNode* root) 

{ 

	return !(root->left) && !(root->right); 
} 

// Creates a min heap of capacity 
// equal to size and inserts all intacter of 
// data[] in min heap. Initially size of 
// min heap is equal to capacity 
struct MinHeap* createAndBuildMinHeap(short int data[],long int freq[], struct HUFFTABLE **address, int size) 
{ 

	struct MinHeap* minHeap = createMinHeap(size); 

	for (int i = 0; i < size; ++i) 
		minHeap->array[i] = newNode(data[i], freq[i], address[i]); 

	// for (int i = 0; i < size; ++i) 
	// 	printf("3h: %d %ld %p\n", minHeap->array[i]->data, minHeap->array[i]->freq, minHeap->array[i]->add);
	
	minHeap->size = size; 
	buildMinHeap(minHeap); 

	
	return minHeap; 
} 


void inorder(struct MinHeapNode* root) 
{ 

	// Assign 0 to left edge and recur 
	if (root->left) { 
		// arr[top] = 0; 
		// printCodes(root->left, arr, top + 1); 
		inorder(root->left);
	} 
	printf("here: %d %ld %p ", root->data, root->freq, root->add);
	if(root->add)
	printf("%d",(*(root->add)).value);
	printf("\n");
	// Assign 1 to right edge and recur 
	if (root->right) { 

		// arr[top] = 1; 
		// printCodes(root->right, arr, top + 1); 
		inorder(root->right);
	} 

 
} 

// The main function that builds Huffman tree 
struct MinHeapNode* buildHuffmanTree(short int data[], long int freq[], struct HUFFTABLE **address, int size) 

{ 
	struct MinHeapNode *left, *right, *top; 

	// Step 1: Create a min heap of capacity 
	// equal to size. Initially, there are 
	// modes equal to size. 
	struct MinHeap* minHeap = createAndBuildMinHeap(data, freq, address, size); 
	
	// Iterate while size of heap doesn't become 1 
	while (!isSizeOne(minHeap)) { 

		// Step 2: Extract the two minimum 
		// freq items from min heap 
		left = extractMin(minHeap); 
		right = extractMin(minHeap); 

		// Step 3: Create a new internal 
		// node with frequency equal to the 
		// sum of the two nodes frequencies. 
		// Make the two extracted node as 
		// left and right children of this new node. 
		// Add this node to the min heap 
		// '$' is a special value for internal nodes, not used 
		top = newNode('$', left->freq + right->freq, NULL); 

		top->left = left; 
		top->right = right; 

		insertMinHeap(minHeap, top); 
	} 

	// inorder(extractMin(minHeap));
	// Step 4: The remaining node is the 
	// root node and the tree is complete. 
	return extractMin(minHeap); 
} 

void printArr(int arr[], int n) 
{ 
	int i; 
	for (i = 0; i < n; ++i) 
		printf("%d", arr[i]); 

	printf("\n"); 

} 

// Prints huffman codes from the root of Huffman Tree. 
// It uses arr[] to store codes 

// A utility function to print an array of size n 

void printCodes(struct MinHeapNode* root, int arr[], int top) 
{ 

	// Assign 0 to left edge and recur 
	if (root->left) { 

		arr[top] = 0; 
		printCodes(root->left, arr, top + 1); 
	} 

	// Assign 1 to right edge and recur 
	if (root->right) { 

		arr[top] = 1; 
		printCodes(root->right, arr, top + 1); 
	} 

	// If this is a leaf node, then 
	// it contains one of the input 
	// intacters, print the intacter 
	// and its code from arr[] 
	if (isLeaf(root)) { 

		// printf("%c: ", root->data); 
		struct HUFFTABLE *p = root->add;
		
		unsigned char *temp = (unsigned char *)malloc(sizeof(unsigned char)*top);
		for(int i=0;i<top;++i){
			temp[i] = arr[i] + '0';
		}
		(*p).bitstring = temp ;
		// printf("h: %d %s\n",(*p).value,(*p).bitstring);	
		// printArr(arr, top); 
	} 
} 

// The main function that builds a 
// Huffman Tree and print codes by traversing 
// the built Huffman Tree 
void HuffmanCodes(short int data[],long int freq[],struct HUFFTABLE **address, int size) 

{ 
	// Construct Huffman Tree 
	struct MinHeapNode* root = buildHuffmanTree(data, freq, address, size); 
	// inorder(root);
	// Print Huffman codes using 
	// the Huffman tree built above 
	int arr[MAX_TREE_HT], top = 0; 

	printCodes(root, arr, top); 
}
 
//
//// Driver program to test above functions 
//int main() 
//{ 
//
//	int arr[] = { 1, 2, 3, 4, 5, 6 }; 
//	int freq[] = { 5, 9, 12, 13, 16, 45 }; 
//
//	int size = sizeof(arr) / sizeof(arr[0]); 
//
//	HuffmanCodes(arr, freq, size); 
//
//	return 0; 
//} 


void form_table(struct HUFFTABLE *table, int num)
{
//	num = 5;
//	
//	table = (struct HUFFTABLE *)malloc(sizeof(struct HUFFTABLE)*6);
//	table[0].value = 123 ; table[0].count = 321212 ;
//	table[1].value = 456 ; table[1].count = 123458 ;	
//	table[2].value = 789 ; table[2].count = 1234847 ;
//	table[3].value = 321 ; table[3].count = 789456 ;
//	table[4].value = 345 ; table[4].count = 4578419 ;
//	table[5].value = 265 ; table[5].count = 4547896 ;
//	
    short int *values = (short int *)malloc(sizeof(short int)*(num+1));
    
	long int  *counts = (long int *)malloc(sizeof(long int)*(num+1));
    
	struct HUFFTABLE **address = (struct HUFFTABLE **)malloc(sizeof(struct HUFFTABLE *) * (num+1));
    
    for(int i = 0;i<=num;++i){
    	values[i] = table[i].value;
		counts[i] = table[i].count;
		address[i] = &table[i];
	}
	
	// for(int i = 0;i<=num;++i){
	//    	printf("%d %ld %s\n", (*address[i]).value, (*address[i]).count, address[i]->bitstring );
	// }
	
	HuffmanCodes(values, counts, address, num+1); 
	
	// for(int i = 0;i<=num;++i){
	//    	printf("%d %ld %s\n", table[i].value,table[i].count, table[i].bitstring );
	// }
	
	
	
    
	
	
   /* Please note that the array named table actually has a
      length of num+1. You can access its elements with
      indices from 0 to num. */
	
	

   /****************************************************************************
   * If the Huffman codes were found, this code could print them out to the
   * screen. I have included it to help you understand what the contents of
   * the data structures in the array table should hold. 
   ****************************************************************************/
   
   /*
   if(VERBOSE){
      for(i=0;i<=num;i++){
         if(table[i].count != 0){
            printf("TABLE[%3d] value=%6d count=%6d bits=%4d bitstring=", i,
               table[i].value, table[i].count, table[i].bits);
            for(p=0;p<table[i].bits;p++) printf("%d", (int)table[i].bitstring[p]);
            printf("\n");
         }
      }
   }
   */
}

//
//void form_table(struct HUFFTABLE *table, int num)
//{
//	num = 5;
//	
//	table = (struct HUFFTABLE *)malloc(sizeof(struct HUFFTABLE)*6);
//	table[0].value = 123 ; table[0].count = 321212 ;
//	table[1].value = 456 ; table[1].count = 123458 ;	
//	table[2].value = 789 ; table[2].count = 1234847 ;
//	table[3].value = 321 ; table[3].count = 789456 ;
//	table[4].value = 345 ; table[4].count = 4578419 ;
//	table[5].value = 265 ; table[5].count = 4547896 ;
//	
//    short int *values = (short int *)malloc(sizeof(short int)*(num+1));
//    
//	long int  *counts = (long int *)malloc(sizeof(long int)*(num+1));
//    
//	struct HUFFTABLE **address = (struct HUFFTABLE **)malloc(sizeof(struct HUFFTABLE *) * (num+1));
//    
//    for(int i = 0;i<=num;++i){
//    	values[i] = table[i].value;
//		counts[i] = table[i].count;
//		address[i] = &table[i];
//	}
//	
//
//	for(int i = 0;i<=num;++i){
//    	printf("%d %ld\n",(*address[i]).value,(*address[i]).count);
//	}
//	
//		
//	int arr[] = {1,0,1};
//	int top = 3;
//	struct HUFFTABLE *p = address[0];
//	
//	(*p).bitstring = (unsigned char *)malloc(sizeof(unsigned char)*top);
//	for(int i=0;i<top;++i){
//		(*p).bitstring[i] = arr[i] + '0';
//	}
//	
//	printf("%s",(*address[0]).bitstring);
//	
//	
//	
//	
//	
//    
//	
//	
//   /* Please note that the array named table actually has a
//      length of num+1. You can access its elements with
//      indices from 0 to num. */
//	
//	
//
//   /****************************************************************************
//   * If the Huffman codes were found, this code could print them out to the
//   * screen. I have included it to help you understand what the contents of
//   * the data structures in the array table should hold. 
//   ****************************************************************************/
//   
//   /*
//   if(VERBOSE){
//      for(i=0;i<=num;i++){
//         if(table[i].count != 0){
//            printf("TABLE[%3d] value=%6d count=%6d bits=%4d bitstring=", i,
//               table[i].value, table[i].count, table[i].bits);
//            for(p=0;p<table[i].bits;p++) printf("%d", (int)table[i].bitstring[p]);
//            printf("\n");
//         }
//      }
//   }
//   */
//}
//
//int main(){
//	form_table(NULL, 0);
//}
