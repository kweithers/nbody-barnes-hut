#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#define SOFTENING 1e-9f
#define THETA 0.6f //9999999.f
#define dt .01f

typedef struct {
	float x,y,vx,vy,Fx,Fy;
} Particle;

typedef struct qtnode_{
	Particle *particle;
	int which_child;
	float size;
	float cm_x,cm_y;
	float total_mass;
	float lb,rb,db,ub;
	struct qtnode_* parent;
	struct qtnode_* child[4];
} QTNode;

void lattice_init(Particle *p, int n)
{
	int i, j;
	float dx, dy;
	int m = (int)sqrt(n);
	dx = dy = 100. / m;
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			p[j + m * i].x = dx * (j + .5);
			p[j + m * i].y = dy * (i + .5);
			p[j + m * i].vx = 0;
			p[j + m * i].vy = 0;
			p[j + m * i].Fx = 0;
			p[j + m * i].Fy = 0;

		}
	}
}

QTNode *create_node(QTNode *parent, int child_index)
{
	QTNode *node = malloc(sizeof(QTNode));
	node->which_child = child_index;
	for (int i = 0; i <= 3; ++i)
		node->child[i] = NULL;
	node->particle = NULL;
	node->parent = parent;
	node->total_mass = 0;
	node->cm_x = 0;
	node->cm_y = 0;

	if (parent == NULL)
	{
		node->size = 100.;
		node->lb = 0.;
		node->rb = 100.;
		node->db = 0.;
		node->ub = 100.;
	}
	else
	{
		node->size = .5 * node->parent->size;
		if (node->which_child == 0)
		{
			node->lb = node->parent->lb;
			node->rb = node->lb + node->size;
			node->ub = node->parent->ub;
			node->db = node->ub - node->size;
		}
		else if (node->which_child == 1)
		{
			node->rb = node->parent->rb;
			node->lb = node->rb - node->size;
			node->ub = node->parent->ub;
			node->db = node->ub - node->size;
		}
		else if (node->which_child == 2)
		{
			node->lb = node->parent->lb;
			node->rb = node->lb + node->size;
			node->db = node->parent->db;
			node->ub = node->db + node->size;
		}
		else if (node->which_child == 3)
		{
			node->rb = node->parent->rb;
			node->lb = node->rb - node->size;
			node->db = node->parent->db;
			node->ub = node->db + node->size;
		}
	}

	return node;
}

QTNode *create_root(QTNode *parent, int child_index, float xmin, float xmax, float ymin, float ymax)
{
	QTNode *node = malloc(sizeof(QTNode));
	node->which_child = child_index;
	for (int i = 0; i <= 3; ++i)
		node->child[i] = NULL;
	node->particle = NULL;
	node->parent = parent;
	node->total_mass = 0;
	node->cm_x = 0;
	node->cm_y = 0;

	node->lb = floor(xmin);
	node->rb = ceil(xmax);
	node->db = floor(ymin);
	node->ub = ceil(ymax);

	float xrange = node->rb - node->lb;
	float yrange = node->ub - node->db;
	float xmid = node->lb + xrange/2.;
	float ymid = node->db + yrange/2.;

	if (xrange > yrange)
	{
		yrange = xrange;
		node->ub = ymid + yrange / 2.;
		node->lb = ymid - yrange / 2.;
		node->size = node->ub - node->lb;
	}
	else if (yrange > xrange)
	{
		xrange = yrange;
		node->lb = xmid - xrange / 2.;
		node->rb = xmid + xrange / 2.;
		node->size = node->rb - node->lb;
	}
	node->size = node->rb - node->lb;


	//node->size = (node->rb - node->lb) * (node->ub - node->lb);

	return node;
}

int in_box(Particle *p, float lb, float rb, float db, float ub)
{
	return (p->x >= lb && p->x <= rb && p->y >= db && p->y <= ub) ? 1 : 0;
}

QTNode *which_child_contains(QTNode *n, Particle *p)
{

	if (in_box(p, n->child[0]->lb, n->child[0]->rb, n->child[0]->db, n->child[0]->ub))
		return n->child[0];
	else if (in_box(p, n->child[1]->lb, n->child[1]->rb, n->child[1]->db, n->child[1]->ub))
		return n->child[1];
	else if (in_box(p, n->child[2]->lb, n->child[2]->rb, n->child[2]->db, n->child[2]->ub))
		return n->child[2];
	else if (in_box(p, n->child[3]->lb, n->child[3]->rb, n->child[3]->db, n->child[3]->ub))
		return n->child[3];
}

int is_leaf(QTNode *n)
{
	if (n == NULL)
		return 0;
	return (n->child[0] == NULL && n->child[1] == NULL && n->child[2] == NULL && n->child[3] == NULL) ? 1 : 0;
}

void qTree_insert(Particle *p, QTNode *n)
{
	QTNode *c;

	//if n is an internal node; i.e. n has a non NULL child
	if (is_leaf(n) == 0)
	{
		c = which_child_contains(n, p);
		n->cm_x = ((n->cm_x * (n->total_mass / (n->total_mass + 1))) + (p->x * (1 / (n->total_mass + 1))));
		n->cm_y = ((n->cm_y * (n->total_mass / (n->total_mass + 1))) + (p->y * (1 / (n->total_mass + 1))));
		n->total_mass += 1.f;
		qTree_insert(p, c);
	}

	// n is a leaf, if n contains a particle;
	else if (n->particle != NULL)
	{
		for (int i = 0; i <= 3; ++i)
			n->child[i] = create_node(n, i);
		//find which child contains the current particle
		c = which_child_contains(n, n->particle);
		//check the insert makes sense
		//printf("Moving particle (%f,%f) to between (%f,%f) and (%f,%f)\n", n->particle->x, n->particle->y, c->lb, c->rb, c->db, c->ub);

		//store the current particle in the correct child
		c->particle = n->particle;

		//update that child's center of mass and total mass
		c->cm_x = ((c->cm_x * (c->total_mass / (c->total_mass + 1))) + (n->particle->x * (1 / (c->total_mass + 1))));
		c->cm_y = ((c->cm_y * (c->total_mass / (c->total_mass + 1))) + (n->particle->y * (1 / (c->total_mass + 1))));
		c->total_mass += 1.f;

		//find which particle contains the new particle
		c = which_child_contains(n, p);
		//update the current node's center of mass and total mass
		n->cm_x = ((n->cm_x * (n->total_mass / (n->total_mass + 1))) + (p->x * (1 / (n->total_mass + 1))));
		n->cm_y = ((n->cm_y * (n->total_mass / (n->total_mass + 1))) + (p->y * (1 / (n->total_mass + 1))));
		n->total_mass += 1.f;

		qTree_insert(p, c);
	}

	//if n is empty; store particle P in this node.
	else
	{
		n->particle = p;
		//set the center of mass and total mass
		n->cm_x = p->x;
		n->cm_y = p->y;
		n->total_mass += 1.f;
		//check if the insert makes sense
		//printf("Inserting (%f,%f) between (%f,%f) and (%f,%f)\n", n->particle->x, n->particle->y, n->lb, n->rb, n->db, n->ub);
	}
}

//calc force on particle k due to ALL particles in node n
void tree_force(Particle *k, QTNode *n)
{
	//distance from particle k to CM of particles in n
	float r = sqrtf((k->x - n->cm_x) * (k->x - n->cm_x) + (k->y - n->cm_y) * (k->y - n->cm_y));
	float D = n->size; //size of n

	if (is_leaf(n) == 1 && n->particle == NULL)
	{	
		k->Fx += 0;
		k->Fy += 0;	//if n contains just one particle - same as previous hw
	}
	else if (is_leaf(n) == 1 && n->particle != NULL)
	{
		float dx = n->total_mass * (n->cm_x - k->x);
		float dy = n->total_mass * (n->cm_y - k->y);
		float distSqr = dx * dx + dy * dy + SOFTENING;
		float invDist = 1.0f / sqrtf(distSqr);
		float invDist2 = invDist * invDist;

		k->Fx += dx * invDist2;
		k->Fy += dy * invDist2;
	}
	//else if n contains more than one particle
	else
	{
		if (D / r < THETA)
		{ //if far enough away, we can approximate
			//printf("Using approximation!\n");
			//printf("Particle is (%f,%f). Node CM is (%f,%f)\n", k->x, k->y, n->cm_x, n->cm_y);
			float dx = n->total_mass * (n->cm_x - k->x);
			float dy = n->total_mass * (n->cm_y - k->y);
			float distSqr = dx * dx + dy * dy + SOFTENING;
			float invDist = 1.0f / sqrtf(distSqr);
			float invDist2 = invDist * invDist;

			k->Fx += dx * invDist2;
			k->Fy += dy * invDist2;
		}
		else
		{ //if close, need to look inside node
			for (int i = 0; i < 4; ++i)
				tree_force(k, n->child[i]);
		}
	}
}

void free_qtree(QTNode *r)
{
	if (is_leaf(r) ==  1)//if leaf, just free itself
	{
		free(r);
	}
	else //if not leaf, free everything below it, then it self
	{
		for (int i=0;i<4;++i)
		{
			free_qtree(r->child[i]);
		}
		free(r);
	}
}

int main(int argc, char **argv)
{
	float biggest_x = 100., biggest_y = 100., smallest_x = 0., smallest_y = 0.;
	int nParticles = 900;
	float next_x,next_y;
	float *buf = malloc(nParticles * sizeof(Particle));
	Particle *p = (Particle *)buf;
	QTNode *root;
	FILE *datafile    = NULL;      /* output file for particle positions */
  	datafile          = fopen("particles.dat","w");
  	fprintf(datafile,"%d %d %d\n", nParticles, 500, 0);

	lattice_init(p, nParticles);
	//smaller number of particles for less clustered visualization
	//p[0].x = .1; p[0].y = .9; p[1].x = .2; p[1].y = .3; p[2].x = .4; p[2].y = .9; p[3].x = .9; p[3].y = .6;
	
	for (int a=0;a<1000;++a)
	{
		root = create_root(NULL, 0, smallest_x,biggest_x,smallest_y,biggest_y);
		for (int i = 0; i < nParticles; ++i)
		{
			qTree_insert(p + i, root);
			//printf("%d DONE\n", i);
		}

		for (int i = 0; i < nParticles; ++i)
		{
			tree_force(p + i, root);
			/* update instantaneous velocity based on force and timestep */
			(p + i)->vx += dt * (p + i)->Fx;
			(p + i)->vy += dt * (p + i)->Fy;
			/*reset Fx and Fy to zero for future loops*/
			(p + i)->Fx = 0;
			(p + i)->Fy = 0;
		}

		/*update positions*/
		for (int i = 0; i < nParticles; ++i)
		{ /* compute new position */
			next_x = (p + i)->x + (p + i)->vx * dt;
			next_y = (p + i)->y + (p + i)->vy * dt;

			//if bounce off right/left wall; reverse x velocity
			if (next_x < 0 || next_x > 100)
			{
				(p+i)->vx *=-1;	
			}	
			//if bounce off top/bottom wall; reverse y velocity
			if (next_y < 0 || next_y > 100)
			{
				(p+i)->vy *=-1;	
			}	
			
			(p + i)->x += (p + i)->vx * dt;
			(p + i)->y += (p + i)->vy * dt;

			if ((p + i)->x > biggest_x)
				biggest_x = (p + i)->x;
			if ((p + i)->x < smallest_x)
				smallest_x = (p + i)->x;
			if ((p + i)->y > biggest_y)
				biggest_y = (p + i)->y;
			if ((p + i)->y < smallest_y)
				smallest_y = (p + i)->y;
		}
		free_qtree(root);
		for (int i = 0;i < nParticles; ++i)
      		fprintf(datafile, "%d %f %f \n",a, p[i].x, p[i].y);
	}
	free(buf);
	fclose(datafile);
	return 0;
}