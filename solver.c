#include <stdio.h>
#include"mmio.h"
#include"mmio.c"
#include<time.h>
FILE *f, *ofp;
long long int memcount = 0;


// A structure to represent an adjacency list node
struct AdjacencyListNode
{
	int dest;
	struct AdjacencyListNode* next;
	double weight;
};

struct checkListNode
{
	int i;
	int j;
	struct checkListNode* next;
};

struct checkListNode* newcheckListNode(int i, int j)
{
	struct checkListNode* newNode =
		(struct checkListNode*) malloc(sizeof(struct checkListNode));
	newNode->i = i;
	newNode->j = j;
	newNode->next = NULL;
	return newNode;
}

// A structure to represent an adjacency liat

struct AdjList
{
	struct AdjacencyListNode* head;  // pointer to head node of list
	struct LexpCellNode* cell;
};

struct Graph
{
	int V;
	struct AdjList* array;
};

struct AdjacencyListNode* newAdjacencyListNode(int dest)
{
	struct AdjacencyListNode* newNode =
		(struct AdjacencyListNode*) malloc(sizeof(struct AdjacencyListNode));
	newNode->dest = dest;
	newNode->next = NULL;
	return newNode;
}

struct Graph* CREATEGRAPH(int V)
{
	struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
	graph->V = V;

	graph->array = (struct AdjList*) calloc(V + 1, sizeof(struct AdjList));

	int i;
	for (i = 1; i <= V; i++)
		graph->array[i].head = NULL;

	return graph;
}

void EDGEADD(struct Graph* graph, int src, int dest, double weight)
{

	if (src != dest)
	{
		struct AdjacencyListNode* newNode = newAdjacencyListNode(dest);
		newNode->next = graph->array[src].head;
		newNode->weight = weight;
		graph->array[src].head = newNode;

		newNode = newAdjacencyListNode(src);
		newNode->weight = weight;
		newNode->next = graph->array[dest].head;
		graph->array[dest].head = newNode;
	}
}

void EDGEADD1(struct Graph* graph, int src, int dest, int weight)
{

	if (src != dest)
	{
		struct AdjacencyListNode* newNode = newAdjacencyListNode(dest);
		newNode->next = graph->array[src].head;
		newNode->weight = weight;
		graph->array[src].head = newNode;

		/*newNode = newAdjacencyListNode(src);
		newNode->weight = weight;
		newNode->next = graph->array[dest].head;
		graph->array[dest].head = newNode;*/
	}
}

void EDGEADD2(struct Graph* graph, int src, int dest, double weight)
{

		struct AdjacencyListNode* newNode = newAdjacencyListNode(dest);
		newNode->next = graph->array[src].head;
		newNode->weight = weight;
		graph->array[src].head = newNode;

		if (src != dest)
		{
			newNode = newAdjacencyListNode(src);
			newNode->weight = weight;
			newNode->next = graph->array[dest].head;
			graph->array[dest].head = newNode;
		}
}

void printGraph(struct Graph* graph)
{
	int v, d;
	for (v = 1; v <= graph->V; v++)
	{
		struct AdjacencyListNode* pCrawl = graph->array[v].head;
		printf("\n Adjacency list of vertex %d\n head ", v);
		while (pCrawl != 0)
		{

			printf("-> %d (wt=%lf)", pCrawl->dest, pCrawl->weight);
			pCrawl = pCrawl->next;
		}
		printf("\n");
		free(pCrawl);
	}
}

struct LexpHeadNode
{
	int flag;
	struct LexpHeadNode* head;
	struct LexpCellNode* next;  // pointer to head node of list
	struct LexpHeadNode* back;
};

struct LexpCellNode
{
	int head;
	struct LexpCellNode* next;
	struct LexpCellNode* back;
	struct LexpHeadNode* flag;
};
struct LexpCellNode* newAdjacencyListNode1(int vertex) //diff2
{
	struct LexpCellNode* newNode =
		(struct LexpCellNode*) calloc(1, sizeof(struct LexpCellNode));
	newNode->head = vertex;
	newNode->next = NULL;
	newNode->back = NULL;
	newNode->flag = NULL;
	return newNode;
}

struct LexpHeadNode* newAdjacencyListNode2()
{
	struct LexpHeadNode* newNode =
		(struct LexpHeadNode*) calloc(1, sizeof(struct LexpHeadNode));
	newNode->head = NULL;
	newNode->next = NULL;
	newNode->back = NULL;
	newNode->flag = 0;
	return newNode;
}

void Lexp(struct Graph* graph, int *alpha, int *alphainv)
{
	int i;
	int y = graph->V;
	struct AdjacencyListNode *w;
	struct LexpHeadNode* RootHeadNode = newAdjacencyListNode2();
	struct LexpHeadNode* FirstHeadNode = newAdjacencyListNode2();
	RootHeadNode->back = NULL;
	RootHeadNode->next = NULL;
	RootHeadNode->head = FirstHeadNode;
	RootHeadNode->flag = 0;
	FirstHeadNode->flag = 0;
	FirstHeadNode->head = NULL;
	FirstHeadNode->back = RootHeadNode;
	struct LexpCellNode* nextptr;
	struct LexpHeadNode* freevar;

	//Making cell nodes and putting adjacent to the first head node.
	for (i = graph->V; i > 0; i--)//diff1
	{
		if (i == graph->V)
		{
			struct LexpCellNode* NewCellNode = newAdjacencyListNode1(i);
			NewCellNode->flag = FirstHeadNode;
			FirstHeadNode->next = NewCellNode;
			graph->array[i].cell = NewCellNode;
		}
		else
		{
			struct LexpCellNode* NewCellNode = newAdjacencyListNode1(i);
			NewCellNode->flag = FirstHeadNode;
			graph->array[i].cell = NewCellNode;
			NewCellNode->back = graph->array[i + 1].cell;
			nextptr = graph->array[i + 1].cell;
			nextptr->next = NewCellNode;
		}
	}


	for (i = graph->V; i > 0; i--)
	{
		// skipping the vertices
		FirstHeadNode = RootHeadNode->head;
		while (FirstHeadNode->next == NULL)
		{
			RootHeadNode->head = FirstHeadNode->head;
			FirstHeadNode->head->back = RootHeadNode;
			freevar = FirstHeadNode;
			FirstHeadNode = FirstHeadNode->head;
			free(freevar); //freeing variables.		
		}

		//Select a new vertex p to delete from the doubly data structure.

		struct LexpCellNode* p;
		struct LexpCellNode* q;
		struct LexpHeadNode* headptr;

		p = FirstHeadNode->next;

		if (p->next != NULL)
		{
			FirstHeadNode->next = p->next;
			p->next->back = NULL;//FirstHeadNode
		}
		else
		{
			FirstHeadNode->next = NULL;
		}

		alpha[i] = p->head;
		alphainv[p->head] = i;

		//update2
		for (w = graph->array[p->head].head; w != NULL; w = w->next)
		{
			if (alphainv[w->dest] == 0)  // check if the order has been assigned.
			{   // delete w from the list .

				q = graph->array[w->dest].cell;

				if (q->back == NULL)
				{
					q->flag->next = q->next;
				}
				else
				{
					q->back->next = q->next;
				}

				if (q->next != NULL)
				{
					q->next->back = q->back;
				}

				headptr = q->flag->back;

				if (headptr->flag == 0)
				{
					struct LexpHeadNode* NewHeadNode = newAdjacencyListNode2();
					NewHeadNode->head = headptr->head;
					headptr->head = NewHeadNode;
					NewHeadNode->head->back = NewHeadNode;
					NewHeadNode->back = headptr;
					NewHeadNode->flag = 1;
					NewHeadNode->next = 0;
					headptr = NewHeadNode;//??
					// NewHeadNode = NewHeadNode->head;//??
				}

				// add sselected cell to the new set.

				q->next = headptr->next;
				q->flag = headptr;
				q->back = NULL;


				if (headptr->next != NULL)
				{
					headptr->next->back = q;

				}
				headptr->next = q;
			}
		}

		free(p);

		headptr = RootHeadNode;

		while (headptr->head != NULL)
		{
			headptr->flag = 0;
			headptr = headptr->head;
		}
	}
	fprintf(ofp, "\nLexP:\n");
	fprintf(ofp, "Alpha Ordering:\n");
	for (i = 1; i <= graph->V; i++)
	{
		fprintf(ofp, "%d->", alpha[i]);
		//fprintf("alphainv: %d\n", alphainv[i]);
	}
	free(RootHeadNode);
	free(FirstHeadNode);
}


struct checkListNode* head1;
struct checkListNode* temp1;
struct checkListNode* head2;
struct checkListNode* temp2;
struct checkListNode* pivot1;
struct checkListNode* pivot2;
int fillin(struct Graph* graph, int *alpha, int *alphainv)
{
	int m = 0;
	int x, j;
	int i = 0;
	int count = 0;
	struct AdjacencyListNode *w;
	struct AdjacencyListNode *u;
	struct AdjacencyListNode *z;
	int *test = (int *)calloc((graph->V) + 1, sizeof(int));
	int *y = (int *)calloc((graph->V) + 1, sizeof(int));
	int *elim = (int *)calloc((graph->V) + 1, sizeof(int));
	struct checkListNode* pivot = (struct checkListNode*) malloc(sizeof(struct checkListNode));
	pivot->i = 0;
	pivot->j = 0;
	pivot->next = NULL;
	pivot1 = pivot;
	/*struct checkListNode* head;
	struct checkListNode* temp;*/

	fprintf(ofp, "\nFill in Edges:\n");

	for (i = 1; i < graph->V; i++)
	{
		int k = graph->V;
		int v = alpha[i];
		//eliminate duplicates in Adj[v] and compute m[v]
		struct AdjacencyListNode *prev = graph->array[v].head;
		for (w = graph->array[v].head; w != NULL; w = w->next)
		{
			if (elim[w->dest] == 0)
			{
				if (test[alphainv[w->dest]])
				{
					//delete w from adj[v]
					y[w->dest] = 0;
					count--;
					prev->next = w->next;
					//graph->array[w->dest].head = graph->array[w->dest].head->next;
					//while ()
				}

				else
				{
					/*if ((y[w->dest]==1))
					{
					printf("Fill in edges: %d, %d\n", x, w->dest);
					}*/
					//chk here
					//k = min(k, alphainv[w->dest]);//define min function

					test[alphainv[w->dest]] = 1;

					if (k > alphainv[w->dest])
					{
						k = alphainv[w->dest];
					}
					prev = w;
				}
			}

			else { prev = w; }

		}

		for (j = 1; j <= graph->V; j++)
		{
			if (y[j] != 0)
			{
				//printf("Fill in edges: %d, %d\n", x, y[j]);
				//fprintf(ofp, "Fillin Edges:");

				if (x < y[j])
				{
					//fprintf(ofp, "%d,%d  ", x, y[j]);
				}
				else
				{
					//fprintf(ofp, "%d,%d  ", y[j], x);
				}

			}
			y[j] = 0;
		}
		m = alpha[k];

		int flag = 0;
		//add reqn fill in edges and reset test
		for (w = graph->array[v].head; w != NULL; w = w->next)
		{
			u = graph->array[m].head;
			//neighb[w->dest] = 1;
			test[alphainv[w->dest]] = 0;
			if (w->dest != m & (elim[w->dest] == 0))
			{
				while (u != NULL)
				{
					if (u->dest == w->dest)
					{
						flag = 1;
					}
					u = u->next;
				}

				if (flag == 0)
				{
					x = m;
					count = count + 1;
					//add w to adj[m]
					EDGEADD1(graph, m, w->dest, 9);
					y[w->dest] = w->dest;
					if (m < w->dest)
					{
						fprintf(ofp, "%d,%d  ", m, w->dest);//print the edge
						head1=newcheckListNode(m, w->dest);
						temp1 = pivot1->next;
						pivot1->next = head1;
						head1->next = temp1;
					}
					else
					{
						fprintf(ofp, "%d,%d  ", w->dest, m);//print the edge
						head1 = newcheckListNode(m, w->dest);
						temp1 = pivot1->next;
						pivot1->next = head1;
						head1->next = temp1;
					}

				}
				else
				{
					flag = 0;
				}
			}
		}
		//printGraph(graph);
		//printf("count is %d\n", count);
		elim[v] = 1;
	}
	printf("count is %d\n", count);
	fprintf(ofp, "\nNumber of Fill-ins: %d\n\n", count);
	free(test);
	free(elim);
	free(y);
	free(graph);
	return count;

}



struct Lexmcell
{
	int dest;
	struct Lexmcell* next;
};

// A structure to represent an adjacency list in the LexM

struct Lexmhead
{
	struct Lexmcell* head;  // pointer of head node in LexM
};



void LexM(struct Graph* graph, int *alpha, int *alphainv)
{

	int i, k, j, p, sel_vertex, m, temp, templ, count;
	float *L = (float *)calloc((graph->V) + 1, sizeof(float));
	int *vmap = (int *)calloc((graph->V) + 1, sizeof(int));
	int *Lint = (int *)calloc((graph->V) + 1, sizeof(int));
	int *vreach = (int *)calloc((graph->V) + 1, sizeof(int));
	struct Lexmhead* reach;
	struct AdjacencyListNode *w;
	struct AdjacencyListNode *z;
	struct Lexmcell* d;


	for (i = 1; i <= graph->V; i++)
	{
		vmap[i] = i;

	}

	for (i = 1; i <= graph->V; i++)
	{
		L[i] = 1;
		alphainv[i] = 0;
	}

	k = 1;

	//loop

	for (i = graph->V; i > 0; i--)
	{
		//select

		sel_vertex = vmap[graph->V];
		vreach[sel_vertex] = 1;
		L[sel_vertex] = 0;
		alphainv[sel_vertex] = i;
		alpha[i] = sel_vertex;

		reach = (struct Lexmhead*)calloc(k + 1, sizeof(struct Lexmhead));

		for (j = 1; j <= k; j++)
		{
			reach[j].head = NULL;
		}

		//all vertices unnumbered

		for (j = 1; j <= graph->V; j++)
		{
			if (alphainv[j] == 0)
			{
				vreach[j] = 0;
			}
		}

		for (w = graph->array[sel_vertex].head; w != NULL; w = w->next)
		{

			// add w to reach(l(w))
			if (alphainv[w->dest] == 0)
			{
				struct Lexmcell* node = (struct Lexmcell*) malloc(sizeof(struct Lexmcell));
				node->dest = w->dest;
				node->next = reach[(int)(L[w->dest])].head;
				reach[(int)(L[w->dest])].head = node;
				vreach[w->dest] = 1;
				L[w->dest] = L[w->dest] + 0.5;
				// Mark {v,w} as an edge of G
			}
		}



		for (j = 1; j <= k; j++)
		{
			while (reach[j].head != NULL)
			{
				//delete a vertex w from reach(j)??
				d = reach[j].head;
				reach[j].head = reach[j].head->next;

				for (z = graph->array[d->dest].head; z != NULL; z = z->next)
				{

					if (vreach[z->dest] != 1)
					{
						vreach[z->dest] = 1;

						if (L[z->dest] > j)
						{
							struct Lexmcell* nodez = (struct Lexmcell*) malloc(sizeof(struct Lexmcell));
							nodez->dest = z->dest;
							nodez->next = reach[(int)(L[z->dest])].head;
							reach[(int)(L[z->dest])].head = nodez;
							L[z->dest] = L[z->dest] + (0.5);
							// mark {v,z} as an edge of G
						}
						else
						{
							struct Lexmcell* nodez = (struct Lexmcell*) malloc(sizeof(struct Lexmcell));
							nodez->dest = z->dest;
							nodez->next = reach[j].head;
							reach[j].head = nodez;
						}
					}
				}
				free(d);
			}
		}

		for (p = 1; p <= graph->V; p++)
		{
			if (L[p] != 0)
			{
				Lint[p] = ((10 * L[p]));
			}
		}


		for (m = 1; m <= graph->V; m++)
		{
			vmap[m] = m;
		}

		free(reach);
		int size = (graph->V) + 1;
		int large = (20 * k) + 1;
		int *semiSorted, *bsort;
		semiSorted = (int*)calloc(size, sizeof(int));
		bsort = (int*)calloc(size, sizeof(int));
		int significantDigit = 1;
		int largestNum = large;

		// Loop until we reach the largest significant digit
		while (largestNum / significantDigit > 0)
		{

			int bucket[10] = { 0 };

			// Counts the number of "keys" or digits that will go into each bucket
			for (p = 1; p < size; p++)
				bucket[(Lint[p] / significantDigit) % 10]++;

			/*
			* Add the count of the previous buckets,
			* Acquires the indexes after the end of each bucket location in the array
			* Works similar to the count sort algorithm
			*/
			for (p = 1; p < 10; p++)
				bucket[p] += bucket[p - 1];

			// Use the bucket to fill a "semiSorted" array
			for (p = size - 1; p > 0; p--)
			{
				int l = --bucket[(Lint[p] / significantDigit) % 10];
				semiSorted[l] = Lint[p];
				bsort[l] = vmap[p];

			}

			for (p = 0; p < size - 1; p++)
				Lint[p + 1] = semiSorted[p];

			for (p = 0; p < size - 1; p++)
				vmap[p + 1] = bsort[p];

			// Move to next significant digit
			significantDigit *= 10;


		}

		m = 1;

		while (Lint[m] == 0)
		{
			m = m + 1;
		}

		temp = 0;
		templ = 0;
		count = 0;

		for (m; m <= graph->V; m++)
		{
			if (Lint[m] == temp)
			{
				Lint[m] = count;
			}
			else
			{
				count = count + 1;
				templ = Lint[m];
				Lint[m] = count;
				temp = templ;
			}

		}

		for (m = 1; m <= graph->V; m++)
		{
			L[vmap[m]] = Lint[m];
		}

		k = Lint[graph->V];
	}

	//printf("\n");
	fprintf(ofp, "\nLexM:\n");
	fprintf(ofp, "Alpha Ordering:\n");
	for (i = 1; i <= graph->V; i++)
	{

		fprintf(ofp, "%d->", alpha[i]);
		//fprintf(ofp, "%d->", alpha1[i]);
	}
	//printf("\n");
	free(vmap);
	free(Lint);
	free(vreach);
	free(L);

}

rowchange(double ** matrix,int* ordermap, int src, int dest,int size)
{
	double c;
	int i;
	for (i = 1; i <= size+1; i++)
	{
		c = matrix[ordermap[src]][i];
		matrix[ordermap[src]][i] = matrix[dest][i];
		matrix[dest][i] = c;

		//printmatrix(matrix, size);

		/*c = matrix[i][ordermap[src]];
		matrix[i][ordermap[src]] = matrix[i][dest];
		matrix[i][dest] = c;
	 */
	}

	//printmatrix(matrix,size);
}

double ** matrixchange(double ** matrix, int* ordermap, int *alpha, int M)
{
	int l1,i,j;
	double tempaks;

	for (i = 1; i<M + 1; i++)
	{
		if (ordermap[i] != alpha[i])
		{
			//l=i;
			for (j = 1; j <= M; j++)
			{
				if (ordermap[j] == alpha[i])
				{
					tempaks = ordermap[j];
					ordermap[j] = ordermap[i];
					ordermap[i] = tempaks;
					break;
				}
			}
			l1 = j;
			tempaks = matrix[i][M+1];
			matrix[i][M + 1] = matrix[l1][M + 1];
			matrix[l1][M + 1] = tempaks;
			for (j = 1; j<M + 1; j++)
			{
				tempaks = matrix[i][j];
				matrix[i][j] = matrix[l1][j];
				matrix[l1][j] = tempaks;
			}
			for (j = 1; j<M + 1; j++)
			{
				tempaks = matrix[j][i];
				matrix[j][i] = matrix[j][l1];
				matrix[j][l1] = tempaks;
			}
		}
	}

	return matrix;
}

colchange(double ** matrix, int* ordermap, int* ordermap1, int src, int dest, int size)
{
	double c;
	int i;
	for (i = 1; i <= size; i++)
	{
		c = matrix[i][ordermap[src]];
		matrix[i][ordermap[src]] = matrix[i][dest];
		matrix[i][dest] = c;

		
	}

	//c = ordermap[src];

	c = ordermap[src];
	ordermap[src] = dest;
	ordermap[ordermap1[dest]] = c;

	c=ordermap1[dest];
	ordermap1[dest]=src;
	ordermap1[src] = c;


	//printmatrix(matrix, size);

	//printf("\n\n");

	/*for (i = 1; i <= size; i++)
	{
		printf("%d ", ordermap[i]);
	}*/
	
}




double * guass(double** array, int size, int* alpha,int* alphainv)
{
	fprintf(ofp, "Inside guass function\n");
	//printmatrix(array, size);

		int i, j, k, n,flag,count,error;
		error = 0;
		count = 0;
		flag = 0;
		float c;
		float sum = 0.0;
		n = size;

		struct checkListNode* pivot = (struct checkListNode*) malloc(sizeof(struct checkListNode));
		pivot->i = 0;
		pivot->j = 0;
		pivot->next = NULL;
		pivot2 = pivot;

		double* x = (double*)calloc(size + 1, sizeof(double));

		//double **matrix1;
		//matrix1 = (double **)calloc(size + 1, sizeof(double*));


		//for (i = 0; i < size + 1; i++)
		//	matrix1[i] = (double*)calloc(size + 2, sizeof(double));

		//for (j = 1; j <= n; j++) /* loop for the generation of upper triangular matrix*/
		//{
		//	for (i = 1; i <= n+1; i++)
		//	{
		//		matrix1[j][i] = array[j][i];
		//	}
		//	
		//}
		////printmatrix(matrix1, size);
		
		for (j = 1; j <= n; j++) /* loop for the generation of upper triangular matrix*/
		{
			for (i = 1; i <= n; i++)
			{
				if (i>j)
				{
					//printf("\n");
					//printmatrix(array, size);

					if (array[j][j] == 0)
					{   
						printf("The GEM stops here as element is zero");
						error = 1;
						goto here;
					}

					c = array[i][j] / array[j][j];
					
					array[i][j] = 0;
					
					for (k = j+1; k <= n + 1; k++)
					{
						if (array[i][k] == 0 & k!=n+1 )
						{
							flag = 1;
						}
						array[i][k] = array[i][k] - c*array[j][k];
						
						if (flag == 1 & k != n + 1  )
						{
							if (array[i][k] !=0 & k>i)
							{
								count++;

								//if (alpha[i] < alpha[k])
								//{
									fprintf(ofp, "(%d,%d %lg) ", alpha[i], alpha[k], array[i][k]);
									head2 = newcheckListNode(alpha[i], alpha[k]);
									temp2 = pivot2->next;
									pivot2->next = head2;
									head2->next = temp2;
					/*			}
								else
								{
									fprintf(ofp, "(%d,%d %lf) ", alpha[k], alpha[i], array[i][k]);
									head2 = newcheckListNode(alpha[k], alpha[i]);
									temp2 = pivot2->next;
									pivot2->next = head2;
									head2->next = temp2;
								}*/
							}

							flag = 0;
						}
					}
				}
			}
			
		}

		//printmatrix(array, size);

		//printmatrix(matrix1, size);
		//count = 0;
		//fprintf(ofp, "\n");
		//for (j = 1; j <= n; j++) /* loop for the generation of upper triangular matrix*/
		//{
		//	for (i = 1; i <= n; i++)
		//	{
		//		if (matrix1[j][i] == 0 & array[j][i] != 0 & i > j)
		//		{
		//			count++;
		//		}

		//	}
		//}

		

		if (array[n][n] != 0)
		{
			x[n] = array[n][n + 1] / array[n][n];
		}
		else
		{
			printf("GEM has stopped here");
		}
		/* this loop is for backward substitution*/
		
		for (i = n - 1; i >= 1; i--)
		{
			sum = 0;
		
			for (j = i + 1; j <= n; j++)
			{
				sum = sum + array[i][j] * x[j];
			}

			if (array[i][i] != 0)
			{
				x[i] = (array[i][n + 1] - sum) / array[i][i];

			}

			else
			{
				printf("GEM has stopped here");
			}
			
		}
	here:
		fprintf(ofp,"\nNon zero element count = %d", count);
		printf("\nNon zero element count = %d \n", count);
		fprintf(ofp,"\nThe solution is: \n");
		if (error == 0)
		{
			for (i = 1; i <= size; i++)
			{
				fprintf(ofp, "x%d=%f\t", i, x[alphainv[i]]); /* x1, x2, x3 are the required solutions*/
			}

			fprintf(ofp, "\n");
		}
		else
		{
			fprintf(ofp, "\nNo Solution exists");
			error = 0;
		}
		int chkcnt=0;
		struct checkListNode* tempchk1;
		struct checkListNode* tempchk2;
		tempchk1 = pivot1;
		tempchk2 = pivot2;
		while (pivot1 != NULL)
		{
			while (pivot2 != NULL)
			{
				if ((pivot1->i == pivot2->i & pivot1->j == pivot2->j )|| (pivot1->i == pivot2->j & pivot1->j == pivot2->i))
				{
					chkcnt++;
					pivot1->i=0;
					pivot1->j=0;
     			}
				pivot2 = pivot2->next;
			}
			pivot2 = tempchk2;
			pivot1 = pivot1->next;
		}
		pivot1 = tempchk1;

		/*printf("\n");
		while (pivot1 != NULL)
		{
			if (pivot1->i != 0 & pivot1->j != 0)
			{
				printf("%d %d ", pivot1->i, pivot1->j);
			}
			pivot1 = pivot1->next;
		}*/

		printf("\n");
		chkcnt--;
		fprintf(ofp, " Matched=%d Notmatched=%d\n", chkcnt, count-chkcnt);	
		//printmatrix(array, size);
 return x;
}

printmatrix(double ** array, int size)
{
	fprintf(ofp,"\n\n");
	int i, j,temp;
	for (i = 1; i <= size; i++)
	{
		for (j = 1; j <= size+1; j++)
		{
			//temp = (int)(array[i][j]);
			fprintf(ofp,"%lg  ", array[i][j]);
		}
		fprintf(ofp,"\n");
	}
}
int main(int argc, char *argv[])
{
	clock_t t1, t2;
	t1 = clock();
	int ret_code;
	char *inFileString, *outFileString;
	MM_typecode matcode;
	inFileString = argv[1];
	outFileString = argv[2];
	//FILE *f, *ofp;
	int M, N, nz;
	int i, *I, *J;
	double *val;
	/*if (argc < 2)
	{
	fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
	exit(1);
	}
	else
	{
	if ((f = fopen(argv[1], "r")) == NULL)
	exit(1);
	}*/

	f = fopen(inFileString, "r");

	if (f == NULL)
	{
		printf("File null\n");
		exit(1);
	}
	if (f != NULL)
	{
		printf("FIle is not null\n");
	}


	if (mm_read_banner(f, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
		mm_is_sparse(matcode))
	{
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	/* find out size of sparse matrix .... */

	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
		exit(1);


	/* reseve memory for matrices */

	I = (int *)malloc(nz * sizeof(int));
	memcount = memcount + (nz * sizeof(int));
	J = (int *)malloc(nz * sizeof(int));
	memcount = memcount + (nz * sizeof(int));
	val = (double *)malloc(nz * sizeof(double));
	memcount = memcount + (nz * sizeof(double));


	/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
	/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
	/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

	for (i = 0; i < nz; i++)
	{
		fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
		//I[i]--;  /* adjust from 1-based to 0-based */
		//J[i]--;
	}

	for (i = 0; i < nz; i++)
	{
		if (val[i] == 0)
		{
			I[i] = 0;
			J[i] = 0;
		}
	}

	int j;
	for (i = 0; i < nz; i++)
	{
		for (j = i + 1; j < nz; j++)
		{
			if ((I[i] == J[j]) & (J[i] == I[j]) & I[i]!=0 & J[i]!=0 )
			{
				if (val[i] >= val[j])
				{
					val[i] = val[j];
				}

				I[j] = 0;
				J[j] = 0;
				val[j] = 0;
			
				
			}

		}
	}

	if (f != stdin) fclose(f);

	/************************/
	/* now write out matrix */
	/************************/

	mm_write_banner(stdout, matcode);
	mm_write_mtx_crd_size(stdout, M, N, nz);
	struct Graph* graph = CREATEGRAPH(M);
	//memcount = memcount + ((N + 1)*sizeof(struct AdjacencyListNode));
	struct Graph* graph2 = CREATEGRAPH(M);
	//memcount = memcount + ((N + 1)*sizeof(struct AdjacencyListNode));
	struct Graph* graph3 = CREATEGRAPH(M);
	//memcount = memcount + ((N + 1)*sizeof(struct AdjacencyListNode));

	
	/*struct Graph* graph4 = CREATEGRAPH(M);
	memcount = memcount + ((N + 1)*sizeof(struct AdjacencyListNode));*/

	for (i = 0; i < nz; i++)
	{
		EDGEADD(graph, I[i], J[i], val[i]);
		EDGEADD(graph2, I[i], J[i], val[i]);
		EDGEADD(graph3, I[i], J[i], val[i]);
		//fprintf(stdout, "%d %d %20.19g\n", I[i] + 1, J[i] + 1, val[i]);
	}
	//fprintf(stdout, "%d %d %20.19g\n", I[i] + 1, J[i] + 1, val[i]);
	//printGraph(graph);

	char outputFilename[] = "major.txt";
	ofp = fopen(outFileString, "w");
	if (ofp == NULL)
	{
		printf("file is not read properly\n");
	}
	fprintf(ofp, "Name of the file %s \n" ,inFileString);
	fprintf(ofp, "Size: %d * %d\nNonZero Elements: %d\n", M, N, nz);

	double* b;
	b = (double *)calloc(M + 1, sizeof(double));
	memcount = memcount + ((M + 1)*sizeof(double));


	fprintf(ofp, "\nB matrix values are \n");

	for (i = 1; i <= M; i++)
	{
		if (mm_is_integer(matcode))
		{
			b[i] = rand() % 10;
			fprintf(ofp, "%d ", b[i]);
		}
		else
		{
			b[i] = rand() % 200;
			b[i] = b[i] / 1.0000;
			b[i] = b[i] / RAND_MAX;
			fprintf(ofp, "%lg ", b[i]);
		}
		}

	fprintf(ofp, "\n");
	

	//////////////////////////////////////////////lexp////////////////////

	int	*alpha1lexp;
	int	*alpha2lexp;
	alpha1lexp = (int *)calloc(M + 1, sizeof(int));
	//memcount = memcount + ((N + 1)*sizeof(int));
	alpha2lexp = (int *)calloc(M + 1, sizeof(int));
	//memcount = memcount + ((N + 1)*sizeof(int));
	int count = 0;
	// ordering of elimination
	//memcount = memcount + (3 * (N + 1)*sizeof(int));
	Lexp(graph, alpha1lexp, alpha2lexp);
	//memcount = memcount + ((N + 1)*sizeof(struct LexpCellNode));
	printf("\n LexP number of Fillin ");
	count = fillin(graph, alpha1lexp, alpha2lexp);

	double **matrix;
	// allocate the 2d array
	matrix = (double **)calloc(M + 1, sizeof(double*));
	for (i = 0; i < M + 1; i++)
		matrix[i] = (double*)calloc(M + 2, sizeof(double));
	memcount = memcount + ((M*M)*sizeof(double));

	//printmatrix(matrix, M);
	for (i = 0; i < nz; i++)
	{
		matrix[I[i]][J[i]] = val[i];
		matrix[J[i]][I[i]] = val[i];
	}


	for (i = 1; i <= M; i++)
	{
		matrix[i][M + 1] = b[i];
	}

	//printmatrix(matrix, M);

	int	*ordermap;


	ordermap = (int *)calloc(M + 1, sizeof(int));
	memcount = memcount + ((M + 1)*sizeof(int));

	for (i = 1; i <= M; i++)
	{
		ordermap[i] = i;
	}

	//do swapping
	
	matrix=matrixchange(matrix, ordermap, alpha1lexp, M);

	double* l;

	//printf("\n");

	l=guass(matrix, M,alpha1lexp,alpha2lexp);
	memcount = memcount + ((M + 1)*sizeof(double));


	
	
	free(matrix);
	free(ordermap);
	//free(ordermap1);
	free(alpha1lexp);
	free(alpha2lexp);
	//free(pivot1);
	//free(pivot2);


///////////////////////////////////////lexm////////////////////////

	
	
	int	*alpha1lexm;
	int	*alpha2lexm;
	alpha1lexm = (int *)calloc(M + 1, sizeof(int));
	//memcount = memcount + ((N + 1)*sizeof(int));
	alpha2lexm = (int *)calloc(M + 1, sizeof(int));
	//memcount = memcount + ((N + 1)*sizeof(int));

	LexM(graph2, alpha1lexm, alpha2lexm);
	//memcount = memcount + ((N + 1)*sizeof(struct Lexmcell));
	printf("\n LexM number of Fillin ");
	count = fillin(graph2, alpha1lexm, alpha2lexm);
	//}

	
	// allocate the 2d array
	matrix = (double **)calloc(M + 1, sizeof(double*));
	for (i = 0; i < M + 1; i++)
		matrix[i] = (double*)calloc(M + 2, sizeof(double));
	//printmatrix(matrix, M);
	for (i = 0; i < nz; i++)
	{
		matrix[I[i]][J[i]] = val[i];
		matrix[J[i]][I[i]] = val[i];
	}

	for (i = 1; i <= M; i++)
	{
		matrix[i][M + 1] = b[i];
	}

	//printmatrix(matrix, M);

	ordermap = (int *)calloc(M + 1, sizeof(int));

	for (i = 1; i <= M; i++)
	{
		ordermap[i] = i;
	}


	matrix=matrixchange(matrix, ordermap, alpha1lexm, M);
	l=guass(matrix, M,alpha1lexm,alpha2lexm);
	free(matrix);
	free(ordermap);
	free(alpha1lexm);
	free(alpha2lexm);
	//free(ordermap1);


	////////////////////////////original////////////////////////////////

	int	*alpha1ori;
	int	*alpha2ori;
	alpha1ori = (int *)calloc(M + 1, sizeof(int));
	//memcount = memcount + ((N + 1)*sizeof(int));
	alpha2ori = (int *)calloc(M + 1, sizeof(int));
	//memcount = memcount + ((N + 1)*sizeof(int));
	fprintf(ofp, "\nOriginal Ordering:\n");
	for (i = 1; i <= graph3->V; i++)
	{
		alpha1ori[i] = i;
		alpha2ori[i] = i;
		fprintf(ofp, "%d->", alpha1ori[i]);
	}
	printf("\n Original ordering number of Fillin ");
	count = fillin(graph3, alpha1ori, alpha2ori);

	// allocate the 2d array
	matrix = (double **)calloc(M + 1, sizeof(double*));
	for (i = 0; i < M + 1; i++)
		matrix[i] = (double*)calloc(M + 2, sizeof(double));
	//printmatrix(matrix, M);
	for (i = 0; i < nz; i++)
	{
		matrix[I[i]][J[i]] = val[i];
		matrix[J[i]][I[i]] = val[i];
	}

	for (i = 1; i <= M; i++)
	{
		matrix[i][M + 1] = b[i];
	}

	//printmatrix(matrix, M);

	ordermap = (int *)calloc(M + 1, sizeof(int));

	for (i = 1; i <= M; i++)
	{
		ordermap[i] = i;
	}

	matrix = matrixchange(matrix, ordermap, alpha1ori, M);
	l=guass(matrix, M,alpha1ori,alpha2ori);
	free(matrix);
	free(ordermap);

	fprintf(ofp, "\n");
	fclose(ofp);
	t2 = clock();
	float diff = (((float)t2 - (float)t1)*0.000001);
	printf("\n \n Total Execution time %f seconds \n", diff);
	printf("\n \n Total Memory used is %llu Kilobytes \n", memcount / 1024);
	return 0;
}