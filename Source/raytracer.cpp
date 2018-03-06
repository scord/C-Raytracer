#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <limits>
#include <glm/gtx/rotate_vector.hpp>
#include "glm/ext.hpp"
#include <fstream>
#include <sstream>

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES

																								         */

mat3 R = glm::mat3(0.1f);

bool printing = false;

int frame = 0;

const float yaw = 0.1;

const int lenseD = 1;

const float lenseSpacing = 0.005;
const int scale = 5;
const int SCREEN_WIDTH =  100;
const int SCREEN_HEIGHT = 100;

bool spatialOn = true;                

const int splitDepth = 5;
const int startingSplit = 0; // 0 = x, 1 = y, 2 = z
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
float focalLength = SCREEN_HEIGHT/1.2;
float focalPoint = 1.3;
float SCREEN_HEIGHT_HALF = SCREEN_HEIGHT/2;
float SCREEN_WIDTH_HALF = SCREEN_WIDTH/2;
vec3 cameraPos(0,0,-2);
vec3 lightPos( 0, -0.5, -0.7 );
vec3 lightColor = 14.0f * vec3( 1, 1, 1 );
vec3 specColor = vec3(1,1,1);
const float pi4 = 4.0f*3.1416f;
int depthS = 0;
int xTest = SCREEN_WIDTH_HALF;//24;
int yTest = SCREEN_HEIGHT_HALF;//44;

int cX;
int cY;

int minSpatialLeafNodeSize = 1;

bool p = false;

struct WorldLimits {
    float maxX;
    float minX;

    float maxY;
    float minY;

    float maxZ;
    float minZ;
};

struct PolyTree {
    PolyTree *LTree;
    PolyTree *RTree;

    vector<int> triangles;
};

WorldLimits *world;
PolyTree *spatialTree;

float nsq;
float n = 1.0f;

float n2sq;
float n2 = 2.0f;

	vec3 ***im;

vec3 indirectLight = 0.1f * vec3( 1, 1, 1 );

struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
	vec3 color;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
bool ClosestIntersection(
	vec3 start,
	vec3 dir,
	const vector<Triangle> &triangles,
	Intersection &closestIntersection );

vec3 DirectLight( const Intersection& i );


struct RayIntersections {
  vec3 points[2];
  float dists[2];
  int hits;
};

void savePos(){
  ofstream myfile;
  myfile.open ("pos.data");


  myfile << cameraPos.x;
  myfile  <<  " ";

  myfile << cameraPos.y;
  myfile  <<  " ";

  myfile << cameraPos.z;
  myfile  <<  " " << endl;

  myfile << lightPos.x;
  myfile  <<  " ";

  myfile << lightPos.y;
  myfile  <<  " ";

  myfile << lightPos.z;
  myfile  <<  " " << endl;

  for (int x = 0; x < 3;x++){
        for (int y = 0; y < 3;y++){
            myfile << R[x][y];
            myfile <<  " ";
        }
  }

  myfile << endl;

  myfile.close();
}

void readPos(){
  ifstream myfile;
  myfile.open ("pos.data");

  if (myfile.is_open()){

      string buf; // Have a buffer string

      string line;

      getline(myfile, line);

      stringstream ss(line); // Insert the string into a stream


      ss >> cameraPos.x;


      ss >> cameraPos.y;



      ss >> cameraPos.z;

      getline(myfile, line);

      ss.clear();
	  ss.str(line); // Insert the string into a stream


      ss >> lightPos.x;

      ss >> lightPos.y;

      ss >> lightPos.z;

      getline(myfile, line);


      ss.clear();
	  ss.str(line); // Insert the string into a stream

      for (int x = 0; x < 3;x++){
            for (int y = 0; y < 3;y++){
                      ss >> R[x][y];

            }
      }

      myfile.close();

  }
}

bool squareContainsPoint(float pointX, float pointY, float minX, float minY, float maxX, float maxY){
    if (p){
        cout << pointX << " " << minX << " " << maxX << endl;
        cout << pointY << " " << minY << " " << maxY << endl;
    }
    return (pointX >= minX && pointX <= maxX && pointY >= minY && pointY <= maxY);
}

bool boxContainsPoint(vec3 point, WorldLimits *wl){
    return (squareContainsPoint(point.x, point.y, wl->minX, wl->minY, wl->maxX, wl->maxY) &&
                squareContainsPoint(point.z, point.y, wl->minZ, wl->minY, wl->maxZ, wl->maxY) &&
                squareContainsPoint(point.x, point.z, wl->minX, wl->minZ, wl->maxX, wl->maxZ));
}

float sign(float v0x, float v0y, float v1x, float v1y, float v2x, float v2y)
{
    return (v0x - v2x) * (v1y - v2y) - (v1x - v2x) * (v0y - v2y);
}

bool triangleContainsPoint(float pointX, float pointY, float v0x, float v0y, float v1x, float v1y, float v2x, float v2y){

    bool b1, b2, b3;

    b1 = sign(pointX, pointY, v0x, v0y, v1x, v1y) < 0.0f;
    b2 = sign(pointX, pointY, v1x, v1y, v2x, v2y) < 0.0f;
    b3 = sign(pointX, pointY, v2x, v2y, v0x, v0y) < 0.0f;

    return ((b1 == b2) && (b2 == b3));
}

bool lineIntersectsBox(float startX, float startY, float endX, float endY, float minX, float minY, float maxX, float maxY){
	return !(((startY > maxY) && (endY > maxY)) || ((startY < minY) && (endY < minY)) || ((startX > maxX) && (endX > maxX)) || ((startX < minX) && (endX < minX)));
}

bool triangleBoxIntersect(Triangle tri, WorldLimits wl, int i, bool target){

    bool ansZ = (squareContainsPoint(tri.v0.x, tri.v0.y, wl.minX, wl.minY, wl.maxX, wl.maxY) ||
                    squareContainsPoint(tri.v1.x, tri.v1.y, wl.minX, wl.minY, wl.maxX, wl.maxY) ||
                    squareContainsPoint(tri.v2.x, tri.v2.y, wl.minX, wl.minY, wl.maxX, wl.maxY));

	if (ansZ == false){
			ansZ = (lineIntersectsBox(tri.v0.x, tri.v0.y, tri.v1.x, tri.v1.y, wl.minX, wl.minY, wl.maxX, wl.maxY) ||
					lineIntersectsBox(tri.v1.x, tri.v1.y, tri.v2.x, tri.v2.y, wl.minX, wl.minY, wl.maxX, wl.maxY) ||
					lineIntersectsBox(tri.v2.x, tri.v2.y, tri.v0.x, tri.v0.y, wl.minX, wl.minY, wl.maxX, wl.maxY));
	} else if (ansZ == false){
        ansZ = (triangleContainsPoint(wl.minX, wl.minY, tri.v0.x, tri.v0.y, tri.v1.x, tri.v1.y, tri.v2.x, tri.v2.y) ||
                    triangleContainsPoint(wl.maxX, wl.minY, tri.v0.x, tri.v0.y, tri.v1.x, tri.v1.y, tri.v2.x, tri.v2.y) ||
                    triangleContainsPoint(wl.minX, wl.maxY, tri.v0.x, tri.v0.y, tri.v1.x, tri.v1.y, tri.v2.x, tri.v2.y) ||
                    triangleContainsPoint(wl.maxX, wl.maxY, tri.v0.x, tri.v0.y, tri.v1.x, tri.v1.y, tri.v2.x, tri.v2.y));
    }

    bool ansY = (squareContainsPoint(tri.v0.x, tri.v0.z, wl.minX, wl.minZ, wl.maxX, wl.maxZ) ||
                    squareContainsPoint(tri.v1.x, tri.v1.z, wl.minX, wl.minZ, wl.maxX, wl.maxZ) ||
                    squareContainsPoint(tri.v2.x, tri.v2.z, wl.minX, wl.minZ, wl.maxX, wl.maxZ));

    if (ansY == false){
			ansY = (lineIntersectsBox(tri.v0.x, tri.v0.z, tri.v1.x, tri.v1.z, wl.minX, wl.minZ, wl.maxX, wl.maxZ) ||
					lineIntersectsBox(tri.v1.x, tri.v1.z, tri.v2.x, tri.v2.z, wl.minX, wl.minZ, wl.maxX, wl.maxZ) ||
					lineIntersectsBox(tri.v2.x, tri.v2.z, tri.v0.x, tri.v0.z, wl.minX, wl.minZ, wl.maxX, wl.maxZ));
	} else if (ansY == false){
        ansY = (triangleContainsPoint(wl.minX, wl.minZ, tri.v0.x, tri.v0.z, tri.v1.x, tri.v1.z, tri.v2.x, tri.v2.z) ||
                    triangleContainsPoint(wl.maxX, wl.minZ, tri.v0.x, tri.v0.z, tri.v1.x, tri.v1.z, tri.v2.x, tri.v2.z) ||
                    triangleContainsPoint(wl.minX, wl.maxZ, tri.v0.x, tri.v0.z, tri.v1.x, tri.v1.z, tri.v2.x, tri.v2.z) ||
                    triangleContainsPoint(wl.maxX, wl.maxZ, tri.v0.x, tri.v0.z, tri.v1.x, tri.v1.z, tri.v2.x, tri.v2.z));
    }

    bool ansX = (squareContainsPoint(tri.v0.y, tri.v0.z, wl.minY, wl.minZ, wl.maxY, wl.maxZ) ||
                    squareContainsPoint(tri.v1.y, tri.v1.z, wl.minY, wl.minZ, wl.maxY, wl.maxZ) ||
                    squareContainsPoint(tri.v2.y, tri.v2.z, wl.minY, wl.minZ, wl.maxY, wl.maxZ));

    if (ansX == false){
			ansX = (lineIntersectsBox(tri.v0.z, tri.v0.y, tri.v1.z, tri.v1.y, wl.minZ, wl.minY, wl.maxZ, wl.maxY) ||
					lineIntersectsBox(tri.v1.z, tri.v1.y, tri.v2.z, tri.v2.y, wl.minZ, wl.minY, wl.maxZ, wl.maxY) ||
					lineIntersectsBox(tri.v2.z, tri.v2.y, tri.v0.z, tri.v0.y, wl.minZ, wl.minY, wl.maxZ, wl.maxY));
	} else if (ansX == false){
        ansX = (triangleContainsPoint(wl.minY, wl.minZ, tri.v0.y, tri.v0.z, tri.v1.y, tri.v1.z, tri.v2.y, tri.v2.z) ||
                    triangleContainsPoint(wl.maxY, wl.minZ, tri.v0.y, tri.v0.z, tri.v1.y, tri.v1.z, tri.v2.y, tri.v2.z) ||
                    triangleContainsPoint(wl.minY, wl.maxZ, tri.v0.y, tri.v0.z, tri.v1.y, tri.v1.z, tri.v2.y, tri.v2.z) ||
                    triangleContainsPoint(wl.maxY, wl.maxZ, tri.v0.y, tri.v0.z, tri.v1.y, tri.v1.z, tri.v2.y, tri.v2.z));
    }

    if (target){ cout << "ansX " << ansX << " ansY " << ansY << " ansZ " << ansZ << endl; }

    return (ansX && ansY && ansZ);

}


void worldSplit(PolyTree *t, int splitDirection, int depth, WorldLimits *wl, vector<bool> path){



    if (depth == 0){
        for (int i = 0; i<t->triangles.size();i++){
            cout << "i " << t->triangles[i] << endl;
        }
        return;
    }



    float mid;

    t->LTree = new PolyTree();
    t->RTree = new PolyTree();

    WorldLimits leftWorldLimits;

    WorldLimits rightWorldLimits;

    if (splitDirection == 0){
        mid = (wl->maxX + wl->minX)/2;

        leftWorldLimits = {
        mid,
        wl->minX,

        wl->maxY,
        wl->minY,

        wl->maxZ,
        wl->minZ

        };

        rightWorldLimits = {
        wl->maxX,
        mid,

        wl->maxY,
        wl->minY,

        wl->maxZ,
        wl->minZ

        };

    } else if (splitDirection == 1){

        mid = (wl->maxY + wl->minY)/2;

        leftWorldLimits = {
        wl->maxX,
        wl->minX,

        mid,
        wl->minY,

        wl->maxZ,
        wl->minZ

        };

        rightWorldLimits = {
        wl->maxX,
        wl->minX,

        wl->maxY,
        mid,

        wl->maxZ,
        wl->minZ

        };

    } else if (splitDirection == 2){

        mid = (wl->maxZ + wl->minZ)/2;

        leftWorldLimits = {
        wl->maxX,
        wl->minX,

        wl->maxY,
        wl->minY,

        mid,
        wl->minZ

        };

        rightWorldLimits = {
        wl->maxX,
        wl->minX,

        wl->maxY,
        wl->minY,

        wl->maxZ,
        mid

        };

    }
    cout << path.size() << endl;


    for (int i = 0; i < t->triangles.size(); i++ ){

        Triangle tri = triangles[t->triangles[i]];


        if (triangleBoxIntersect(tri, leftWorldLimits, t->triangles[i], false)){
            t->LTree->triangles.push_back(t->triangles[i]);
        }

    }


    for (int i = 0; i < t->triangles.size(); i++ ){

        Triangle tri = triangles[t->triangles[i]];

            bool target = false;

    if (path.size() == 2) { target = (path[0] && path[1] && t->triangles[i] == 21); }

        if (triangleBoxIntersect(tri, rightWorldLimits, t->triangles[i], target)){
            t->RTree->triangles.push_back(t->triangles[i]);
        }

    }

    if (t->LTree->triangles.size() > minSpatialLeafNodeSize){
            cout << " left " << depth << endl;
            path.push_back(true);
        worldSplit(t->LTree, (splitDirection+1)%3, depth-1, &leftWorldLimits, path);
        path.pop_back();
    }
    if (t->RTree->triangles.size() > minSpatialLeafNodeSize){
            cout << " right " << depth << endl;
            path.push_back(false);
        worldSplit(t->RTree, (splitDirection+1)%3, depth-1, &rightWorldLimits, path);
        path.pop_back();
    }
}

PolyTree* treeMake(){

    PolyTree *t = new PolyTree();

    WorldLimits *wl = new WorldLimits();

    wl->maxX = std::numeric_limits<float>::min();
    wl->minX = std::numeric_limits<float>::max();

    wl->maxY = std::numeric_limits<float>::min();
    wl->minY = std::numeric_limits<float>::max();

    wl->maxZ = std::numeric_limits<float>::min();
    wl->minZ = std::numeric_limits<float>::max();



    for (int i = 0; i < triangles.size(); i++ )
	{
        t->triangles.push_back(i);

        // Xs

        if (wl->maxX < triangles[i].v0.x){
            wl->maxX = triangles[i].v0.x;
        }
        if (wl->maxX < triangles[i].v1.x){
            wl->maxX = triangles[i].v1.x;
        }
        if (wl->maxX < triangles[i].v2.x){
            wl->maxX = triangles[i].v2.x;
        }
        if (wl->minX > triangles[i].v0.x){
            wl->minX = triangles[i].v0.x;
        }
        if (wl->minX > triangles[i].v1.x){
            wl->minX = triangles[i].v1.x;
        }
        if (wl->minX > triangles[i].v2.x){
            wl->maxX = triangles[i].v2.x;
        }

        // Ys

        if (wl->maxY < triangles[i].v0.y){
            wl->maxY = triangles[i].v0.y;
        }
        if (wl->maxY < triangles[i].v1.y){
            wl->maxY = triangles[i].v1.y;
        }
        if (wl->maxY < triangles[i].v2.y){
            wl->maxY = triangles[i].v2.y;
        }
        if (wl->minY > triangles[i].v0.y){
            wl->minY = triangles[i].v0.y;
        }
        if (wl->minY > triangles[i].v1.y){
            wl->minY = triangles[i].v1.y;
        }
        if (wl->minY > triangles[i].v2.y){
            wl->maxY = triangles[i].v2.y;
        }

        // Zs

        if (wl->maxZ < triangles[i].v0.z){
            wl->maxZ = triangles[i].v0.z;
        }
        if (wl->maxZ < triangles[i].v1.z){
            wl->maxY = triangles[i].v1.z;
        }
        if (wl->maxZ < triangles[i].v2.z){
            wl->maxZ = triangles[i].v2.z;
        }
        if (wl->minZ > triangles[i].v0.z){
            wl->minZ = triangles[i].v0.z;
        }
        if (wl->minZ > triangles[i].v1.z){
            wl->minZ = triangles[i].v1.z;
        }
        if (wl->minZ > triangles[i].v2.z){
            wl->maxZ = triangles[i].v2.z;
        }

    }

    wl->maxX += 0.0001f;
    wl->minX -= 0.0001f;
    wl->maxY += 0.0001f;
    wl->minY -= 0.0001f;
    wl->maxZ += 0.0001f;
    wl->minZ -= 0.0001f;

    world = wl;

    vector<bool> path;

    worldSplit(t, startingSplit, splitDepth, wl, path);

    return t;

}


bool contains( vector<int> tris, int j){

    for (int i = 0;i<tris.size();i++){
         if (tris[i] == j){ return true;}

    }
    return false;
}

void traverse(PolyTree* t, vector<int>* tris){
    if (t->LTree == NULL){

        for (int i = 0;i<t->triangles.size();i++){
            if (!contains(*tris, (t->triangles[i]))){
                tris->push_back(t->triangles[i]);
            }
        }
    } else {
        traverse(t->LTree, tris);
        traverse(t->RTree, tris);
    }
}


int main( int argc, char* argv[] )
{


	screen = InitializeSDL( SCREEN_WIDTH*scale, SCREEN_HEIGHT*scale );
	t = SDL_GetTicks();	// Set start value for timer.
	LoadTestModel ( triangles);
	triangles[21].color = vec3(0.75f, 0.15f, 0.15f );
	nsq = n*n;
	n2sq = n2*n2;

	  im = new vec3**[SCREEN_WIDTH];
  for (int i = 0; i < SCREEN_WIDTH; ++i) {
    im[i] = new vec3*[SCREEN_HEIGHT];

    for (int j = 0; j < SCREEN_HEIGHT; ++j)
      im[i][j] = new vec3[lenseD*lenseD];
  }


    spatialTree = treeMake();

    vector<int> tris;

    //traverse(spatialTree, &tris);

    //cout << tris.size() << endl;

   // for (int i = 0;i<tris.size();i++){
   //    cout << tris[i] << endl;
   // }

    readPos();

	while( NoQuitMessageSDL() )
	{

		Update();
		Draw();
	}

    savePos();

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}



vec3 SoftLight(const Intersection& i)
{

	vec3 averageD = vec3(0,0,0);
	vec3 averageS = vec3(0,0,0);
	vec3 normal = triangles[i.triangleIndex].normal;
	float spec = triangles[i.triangleIndex].spec;
	float refl = triangles[i.triangleIndex].refl;
	float shin = triangles[i.triangleIndex].shin;
	float diff = triangles[i.triangleIndex].diff;
	vec3 v = cameraPos - i.position;
	vec3 color = triangles[i.triangleIndex].color;
	vec3** tex = triangles[i.triangleIndex].texture;
	//cout << sizeof(tex) << endl;

	for (int x = 0; x < n2; x++)
	{
		for (int y = 0; y < n2; y++)
		{
			Intersection closestIntersection;
			Intersection closestReflectIntersection;
			vec3 pos = lightPos;

			pos.x = pos.x + -0.05 + x*0.1f/n;
			pos.y = pos.y + -0.05 + y*0.1f/n;
			float r = glm::length(pos - i.position);
			vec3 direction = glm::normalize(pos - i.position);
			vec3 ref = 2.0f * glm::dot(normal, v) * normal - v;
			bool reflectIntersect = false;
			if (refl > 0)
				reflectIntersect = ClosestIntersection(i.position, ref, triangles, closestReflectIntersection);

			vec3 reflectColor(0,0,0);
			if (reflectIntersect)
				reflectColor = SoftLight(closestReflectIntersection);

			vec3 h = glm::normalize((direction + v)/(glm::sqrt(glm::dot(direction +v,direction +v))));
			vec3 B = (lightColor / (pi4*r*r));


			B = i.color*B;


			float S = glm::pow(glm::dot(normal, h),10);
			bool intersect = ClosestIntersection(i.position, direction, triangles, closestIntersection);
		    float D = glm::dot(normal,direction);

			if (S < 0)
				S = 0;

			if (D < 0)
				D = 0;

			if (intersect && closestIntersection.distance < r)
			{
				D = 0;
				S = 0;
			}

			//cout << spec << endl;
			averageD += diff*B*D;
			averageS += spec*B*S + refl*reflectColor;
		}
	}

    return color*averageD/n2sq + averageS/n2sq + color*indirectLight;
}

vec3 DirectLight(const Intersection& i)
{
	Intersection closestIntersection;


	float r = glm::length(lightPos - i.position);
	vec3 B = lightColor / (4.0f*3.1416f*r*r);
	vec3 normal = glm::normalize(triangles[i.triangleIndex].normal);
	vec3 direction = glm::normalize(lightPos - i.position);

	bool intersect = ClosestIntersection(i.position, direction, triangles, closestIntersection);
    float rn = glm::dot(normal,direction);
    vec3 D = B*rn;

	if (intersect && closestIntersection.distance < r)
		D = vec3(0,0,0);

	vec3 color = triangles[i.triangleIndex].color;
    return color*D + color*indirectLight;


}

bool triangleHitTest(vec3 start, vec3 dir, vector<int> tris, Intersection &closestIntersection, WorldLimits *wl){

    float m = std::numeric_limits<float>::max();

    bool intersectionFound = false;

	for (int j = 0; j < tris.size(); j++ )
	{
		int i = tris[j];
		if (xTest == cX && yTest == cY && printing){
			cout << i << endl;
		}

		vec3 v0 = triangles[i].v0;
		vec3 v1 = triangles[i].v1;
		vec3 v2 = triangles[i].v2;

		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = start - v0;

		mat3 A( -dir, e1, e2 );
		vec3 x = glm::inverse( A ) * b;


				if (i == 21 && xTest == cX && yTest == cY && printing){
						vec3 dest = start + dir*x.x;
			cout << "x.x " << x.x << endl;
			cout << "x.y " << x.y << endl;
			cout << "x.z " << x.z << endl;
						cout << "(" << dest.x;
			cout << " ," << dest.y;
			cout << ", " << dest.z << ")" << endl;
						cout << "X: " << wl->minX << endl;
			cout << wl->maxX << endl;
			cout << "Y: " << wl->minY << endl;
			cout << wl->maxY << endl;
			cout << "Z: " << wl->minZ << endl;
			cout << wl->maxZ << endl;

		}

		if (x.y >= 0 && x.z >= 0 && x.y + x.z <= 1 && x.x >= 0.0001f && x.x < m && (boxContainsPoint(start + dir*x.x, wl)))
		{
		if (xTest == cX && yTest == cY && printing){
			cout << "hit!: " << i << endl;
		}
			intersectionFound = true;

            m = x.x;
            closestIntersection.position = start + dir*x.x;
            closestIntersection.distance = x.x;
            closestIntersection.triangleIndex = i;
			vec3** tex = triangles[i].texture;
			int h = triangles[i].texHeight;

			int w = triangles[i].texWidth;
			closestIntersection.color = triangles[i].texture[(int)(w*x.y)][(int)(h*x.z)];
		}
	}

	return intersectionFound;
}

RayIntersections* rayBoxCollide(vec3 start, vec3 dir, WorldLimits *wl){

    RayIntersections *ri = new RayIntersections();

    int numberOfHits = 0;

    // Zs

    float dZmax = (dir.z != 0) ? (wl->maxZ - start.z)/dir.z : 99;



    glm::vec2 Zmax(start.x + dir.x*dZmax, start.y + dir.y*dZmax);


    if (printing && cX == xTest && cY == yTest){
        cout << dZmax << " " << wl->maxZ << " " << start.z << " " << dir.z << endl;
         cout << "ZMax ri (" << Zmax.x << ", " << Zmax.y << ", " << start.z + dir.z*dZmax<< ")" << endl;
    }

    if (squareContainsPoint(Zmax.x, Zmax.y, wl->minX, wl->minY, wl->maxX, wl->maxY)){
		if (printing && cX == xTest && cY == yTest){
			cout << "hit" << endl;
		}

        ri->points[numberOfHits] = vec3(Zmax.x, Zmax.y, start.z + dir.z*dZmax);
        ri->dists[numberOfHits] = dZmax;
        numberOfHits++;
    }


    float dZmin = (dir.z != 0) ? (wl->minZ - start.z)/dir.z : 99;


    glm::vec2 Zmin(start.x + dir.x*dZmin, start.y + dir.y*dZmin);

    if (printing && cX == xTest && cY == yTest){
        cout << dZmin << " " << wl->minZ << " " << start.z << " " << dir.z << endl;
         cout << "ZMin ri (" << Zmin.x << ", " << Zmin.y << ", " << start.z + dir.z*dZmin<< ")" << endl;
    }

    if (squareContainsPoint(Zmin.x, Zmin.y, wl->minX, wl->minY, wl->maxX, wl->maxY)){
        		if (printing && cX == xTest && cY == yTest){
			cout << "hit" << endl;
		}
        ri->points[numberOfHits] = vec3(Zmin.x, Zmin.y, start.z + dir.z*dZmin);

        ri->dists[numberOfHits] = dZmin;
        numberOfHits++;

        if (numberOfHits == 2){
            ri->hits = numberOfHits;
            return ri;
        }
    }



    // Ys



    float dYmax = (dir.y != 0) ? (wl->maxY - start.y)/dir.y : 99;

    glm::vec2 Ymax(start.x + dir.x*dYmax, start.z + dir.z*dYmax);
                    if (printing && cX == xTest && cY == yTest){
        cout << dYmax << " " << wl->maxY << " " << start.y << " " << dir.y << endl;
         cout << "YMax ri (" << Ymax.x << ", " << start.y + dir.y*dYmax << ", " << Ymax.y << ")" << endl;
    }
    if (squareContainsPoint(Ymax.x, Ymax.y, wl->minX, wl->minZ, wl->maxX, wl->maxZ)){
		if (printing && cX == xTest && cY == yTest){
			cout << "hit" << endl;
		}
        ri->points[numberOfHits] = vec3(Ymax.x, start.y + dir.y*dYmax, Ymax.y);
        ri->dists[numberOfHits] = dYmax;
        numberOfHits++;

        if (numberOfHits == 2){
            ri->hits = numberOfHits;
            return ri;
        }
    }



    float dYmin = (dir.y != 0) ? (wl->minY - start.y)/dir.y : 99;



    glm::vec2 Ymin(start.x + dir.x*dYmin, start.z + dir.z*dYmin);
    if (printing && cX == xTest && cY == yTest){
        cout << dYmin << " " << wl->minY << " " << start.y << " " << dir.y << endl;
         cout << "YMin ri (" << Ymin.x << ", " << start.y + dir.y*dYmin << ", " << Ymin.y << ")" << endl;
    }
    if (squareContainsPoint(Ymin.x, Ymin.y, wl->minX, wl->minZ, wl->maxX, wl->maxZ)){
        		if (printing && cX == xTest && cY == yTest){
			cout << "hit" << endl;
		}
		ri->points[numberOfHits] = vec3(Ymin.x, start.y + dir.y*dYmin, Ymin.y);
        ri->dists[numberOfHits] = dYmin;
        numberOfHits++;

        if (numberOfHits == 2){
            ri->hits = numberOfHits;
            return ri;
        }
    }



    // Xs

    float dXmax = (dir.x != 0) ? (wl->maxX - start.x)/dir.x : 99;


    glm::vec2 Xmax(start.y + dir.y*dXmax, start.z + dir.z*dXmax);

        if (printing && cX == xTest && cY == yTest){
        cout << dXmax << " " << wl->maxX << " " << start.x << " " << dir.x << endl;
         cout << "XMax ri (" << start.x + dir.x*dXmax << ", " << Xmax.x << ", " << Xmax.y << ")" << endl;
    }

    if (squareContainsPoint(Xmax.x, Xmax.y, wl->minY, wl->minZ, wl->maxY, wl->maxZ)){
        		if (printing && cX == xTest && cY == yTest){
			cout << "hit" << endl;
		}
		ri->points[numberOfHits] = vec3(start.x + dir.x*dXmax, Xmax.x, Xmax.y);
        ri->dists[numberOfHits] = dXmax;
        numberOfHits++;

        if (numberOfHits == 2){
            ri->hits = numberOfHits;
            return ri;
        }
    }



    float dXmin = (dir.x != 0) ? (wl->minX - start.x)/dir.x : 99;


    glm::vec2 Xmin(start.y + dir.y*dXmin, start.z + dir.z*dXmin);
        if (printing && cX == xTest && cY == yTest){
        cout << dXmin << " " << wl->minX << " " << start.x << " " << dir.x << endl;
         cout << "XMin ri (" << start.x + dir.x*dXmin << ", " << Xmin.x << ", " << Xmin.y << ")" << endl;
    }

    if (squareContainsPoint(Xmin.x, Xmin.y, wl->minY, wl->minZ, wl->maxY, wl->maxZ)){
        		if (printing && cX == xTest && cY == yTest){
			cout << "hit" << endl;
		}
		ri->points[numberOfHits] = vec3(start.x + dir.x*dXmin, Xmin.x, Xmin.y);
        ri->dists[numberOfHits] = dXmin;
        numberOfHits++;

        if (numberOfHits == 2){
            ri->hits = numberOfHits;
            return ri;
        }
    }



    ri->hits = numberOfHits;
    return ri;

}

bool traverseSpatialTree(PolyTree *tree, WorldLimits *wl, int splitDirection, vec3 start, vec3 dir, Intersection &closestIntersection){
    depthS++;
    //cout << "p " << printing << xTest << " " << cX << endl;

	if (printing && cX == xTest && cY == yTest){
cout << "New call: splitDirection: " << splitDirection << " Depth " << depthS << endl;
	}
    if (tree->LTree == NULL){

        return triangleHitTest(start, dir, tree->triangles, closestIntersection, wl);
    }
bool tmp = printing;

printing = false;

    RayIntersections *r = rayBoxCollide(start, dir, wl);
printing = tmp;
    if (r->hits == 0 || (r->dists[0] < -0.00001f && r->dists[1] < -0.00001f)) {
        return false;
    }

    glm::vec3 point1, point2;



    bool b = (r->dists[0] < r->dists[1]);

    float mid;

    if (b){


        if (r->dists[0] < -0.00001f){
            point1 = start;
            r->dists[0] = 0;

        } else {

            point1 = r->points[0];

        }

        point2 = r->points[1];

    } else {

        if (r->dists[1] < -0.00001f){
            point1 = start;
            r->dists[1] = 0;

        } else {

            point1 = r->points[1];

        }

        point2 = r->points[0];
    }

    if (printing && cX == xTest && cY == yTest){
        cout << "hits: " << r->hits << endl;

        cout << "start (" << start.x << ", " << start.y << ", " << start.z << ")" << endl;
        cout << "point1 (" << point1.x << ", " << point1.y << ", " << point1.z << ")" << endl;
        cout << "dist1: " << r->dists[0] << endl;
        cout << "point2 (" << point2.x << ", " << point2.y << ", " << point2.z << ")" << endl;
        cout << "dist2: " << r->dists[1] << endl << endl << endl;
    }
	float midDelt = 0;
    if (splitDirection == 0){

        mid = (wl->maxX + wl->minX)/2;


        WorldLimits newWl;

        bool left = false, right = false;

        if (point1.x <= mid + midDelt){
            left = true;
        newWl = {
            mid,
            wl->minX,

            wl->maxY,
            wl->minY,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p1 left SD: " << splitDirection << endl;
    }
            if (traverseSpatialTree(tree->LTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection)){
				if (printing && cX == xTest && cY == yTest){
					cout << "return, SD: " << splitDirection << " depth: " << depthS << endl;
				}
                return true;
            }

        }

		if (point1.x >= mid - midDelt){
            right = true;
        newWl = {
            wl->maxX,
            mid,

            wl->maxY,
            wl->minY,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p1 right SD: " << splitDirection << endl;
    }
            if (traverseSpatialTree(tree->RTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection)){
				if (printing && cX == xTest && cY == yTest){
					cout << "return, SD: " << splitDirection << " depth: " << depthS << endl;
				}
                return true;
            }

        }

        if (point2.x <= mid  + midDelt && left == false){

        newWl = {
            mid,
            wl->minX,

            wl->maxY,
            wl->minY,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p2 left SD: " << splitDirection << endl;
    }
            return traverseSpatialTree(tree->LTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection);

        }
		if (point2.x >= mid - midDelt && right == false){

        newWl = {
            wl->maxX,
            mid,

            wl->maxY,
            wl->minY,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p2 right SD: " << splitDirection << endl;
    }
            return traverseSpatialTree(tree->RTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection);

        }

    }

    if (splitDirection == 1){

        mid = (wl->maxY + wl->minY)/2;

        WorldLimits newWl;

        bool left = false, right = false;

        if (point1.y <= mid + midDelt){
            left = true;
        newWl = {
            wl->maxX,
            wl->minX,

            mid,
            wl->minY,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p1 left SD: " << splitDirection << endl;
    }
            if (traverseSpatialTree(tree->LTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection)){
				if (printing && cX == xTest && cY == yTest){
					cout << "return, SD: " << splitDirection << " depth: " << depthS << endl;
				}
                return true;
            }

        }
		if (point1.y >= mid - midDelt){
            right = true;
        newWl = {
            wl->maxX,
            wl->minX,

            wl->maxY,
            mid,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p1 right SD: " << splitDirection << endl;
    }
            if (traverseSpatialTree(tree->RTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection)){
				if (printing && cX == xTest && cY == yTest){
					cout << "return, SD: " << splitDirection << " depth: " << depthS << endl;
				}
                return true;
            }

        }

        if (point2.y <= mid + midDelt && left == false){

        newWl = {
            wl->maxX,
            wl->minX,

            mid,
            wl->minY,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p2 left SD: " << splitDirection << endl;
    }
            return traverseSpatialTree(tree->LTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection);

        }
		if (point2.y >= mid - midDelt && right == false){

        newWl = {
            wl->maxX,
            wl->minX,

            wl->maxY,
            mid,

            wl->maxZ,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p2 right SD: " << splitDirection << endl;
    }
            return traverseSpatialTree(tree->RTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection);

        }

    }

    if (splitDirection == 2){

        mid = (wl->maxZ + wl->minZ)/2;

        WorldLimits newWl;

        bool left = false, right = false;

        if (point1.z <= mid + midDelt){
            left = true;
        newWl = {
            wl->maxX,
            wl->minX,

            wl->maxY,
            wl->minY,

            mid,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p1 left SD: " << splitDirection << endl;
    }
            if (traverseSpatialTree(tree->LTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection)){
				if (printing && cX == xTest && cY == yTest){
					cout << "return, SD: " << splitDirection << " depth: " << depthS << endl;
				}
                return true;
            }

        }
		if (point1.z >= mid - midDelt){
           right = true;
        newWl = {
            wl->maxX,
            wl->minX,

            wl->maxY,
            wl->minY,

            wl->maxZ,
            mid

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p1 right SD: " << splitDirection << endl;
    }
            if (traverseSpatialTree(tree->RTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection)){
				if (printing && cX == xTest && cY == yTest){
					cout << "return, SD: " << splitDirection << " depth: " << depthS << endl;
				}
                return true;
            }

        }

        if (point2.z <= mid + midDelt && left == false){

        newWl = {
            wl->maxX,
            wl->minX,

            wl->maxY,
            wl->minY,

            mid,
            wl->minZ

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p2 left SD: " << splitDirection << endl;
    }
            return traverseSpatialTree(tree->LTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection);

        }
		if (point2.z >= mid - midDelt && right == false){

        newWl = {
            wl->maxX,
            wl->minX,

            wl->maxY,
            wl->minY,

            wl->maxZ,
            mid

        };
    if (printing && cX == xTest && cY == yTest){
    	cout << "p2 right SD: " << splitDirection << endl;
    }
            return traverseSpatialTree(tree->RTree, &newWl, (splitDirection+1)%3, start, dir, closestIntersection);

        }

    }


}


bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, Intersection &closestIntersection)
{
    if (spatialOn){
    depthS = 0;
    //cout << printing << " " << cX << " " << cY << endl;
    	bool ans = traverseSpatialTree(spatialTree, world, startingSplit, start, dir, closestIntersection);

      if (printing && cX == xTest && cY == yTest){
        cout << "done traversing" << endl;
      }

      return ans;
    } else {

	float m = std::numeric_limits<float>::max();

    bool intersectionFound = false;

		if (printing && xTest == cX && yTest == cY){
			cout << "start: " << endl;
		}
	for (int i = 0; i < triangles.size(); i++ )
	{

		vec3 v0 = triangles[i].v0;
		vec3 v1 = triangles[i].v1;
		vec3 v2 = triangles[i].v2;

		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = start - v0;

		mat3 A( -dir, e1, e2 );
		vec3 x = glm::inverse( A ) * b;

		if (x.y >= 0 && x.z >= 0 && x.y + x.z <= 1 && x.x >= 0.0001f && x.x < m)
		{
			if (printing && xTest == cX && yTest == cY){
				cout << "hit! " << i << endl;
			}
			intersectionFound = true;

            m = x.x;
            closestIntersection.position = start + dir*x.x;
            closestIntersection.distance = x.x;
            closestIntersection.triangleIndex = i;
			vec3** tex = triangles[i].texture;
			int h = triangles[i].texHeight;

			int w = triangles[i].texWidth;
			closestIntersection.color = triangles[i].texture[(int)(w*x.y)][(int)(h*x.z)];
		}
	}
	return intersectionFound;

	}
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;

	cout << dt << " (ms)" << endl;

	Uint8* keystate = SDL_GetKeyState( 0 );
	if (keystate[SDLK_w])
	{
		lightPos.z += 0.1;
	}
	if (keystate[SDLK_s])
	{
		lightPos.z -= 0.1;
	}
	if (keystate[SDLK_a])
	{
		lightPos.x -= 0.1;
	}
	if (keystate[SDLK_d])
	{
		lightPos.x += 0.1;
	}
	if (keystate[SDLK_e])
	{
		lightPos.y += 0.1;
	}
	if (keystate[SDLK_q])
	{
		lightPos.y -= 0.1;
	}

	if( keystate[SDLK_f] )
	{
		focalPoint-=0.2;
	// increase focal
	}
	if( keystate[SDLK_g] )
	{
		focalPoint+=0.2;
	// dcrease focal
	}

	if( keystate[SDLK_t] )
	{
		xTest = (xTest+4)%SCREEN_WIDTH;
	// increase focal
	}
	if( keystate[SDLK_y] )
	{
		yTest = (yTest+4)%SCREEN_HEIGHT;
	// dcrease focal
	}

    if( keystate[SDLK_p] )
	{
		spatialOn = !spatialOn;
	}

	    if( keystate[SDLK_r] )
	{
		cameraPos = vec3(0,0,-2);
		lightPos = vec3(0, -0.5, -0.7 );
		R = glm::mat3(0.1f);
	// dcrease focal
	}

	if( keystate[SDLK_UP] )
	{
		cameraPos+= vec3( R[2][0], R[2][1], R[2][2] );
	// Move camera forward
	}
	if( keystate[SDLK_DOWN] )
	{
		cameraPos-= vec3( R[2][0], R[2][1], R[2][2] );
	// Move camera backward
	}
	if( keystate[SDLK_LEFT] )
	{
		R[2] = glm::rotateY(glm::vec3( R[2][0], R[2][1], R[2][2] ), -yaw);
		R[0] = glm::rotateY(glm::vec3( R[0][0], R[0][1], R[0][2] ), -yaw);
	// rotate camera to the left
	}
	if( keystate[SDLK_RIGHT] )
	{
		R[2] = glm::rotateY(glm::vec3( R[2][0], R[2][1], R[2][2] ), yaw);
		R[0] = glm::rotateY(glm::vec3( R[0][0], R[0][1], R[0][2] ), yaw);
	// totate camera to the right
	}
}

void Draw()
{
frame++;
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	vec3 d;
	d.z = focalLength;

	cout << "frame " << frame << " FP " <<  focalPoint << endl;

	for (int a = 0; a < lenseD; a++)
	{
		for (int b = 0; b < lenseD; b++)
		{


	for (int x = 0; x < SCREEN_WIDTH; x++)
	{
		for (int y = x % 2; y < SCREEN_HEIGHT; y+=2)
		{
		cX = x;
		cY = y;
		//cout << cX << " " << cY << endl;
			d.x = x - SCREEN_WIDTH_HALF;
			d.y = y - SCREEN_HEIGHT_HALF;

			vec3 averageC = vec3(0,0,0);
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					vec3 pos = d;

					//cout << "x " << x << " y " << y << " a " << a << " b " << b << endl;
					pos.x = pos.x + -0.5 + i*0.5f/n;
					pos.y = pos.y + -0.5 + j*0.5f/n;
					Intersection closestIntersection;
					if (printing && cX == xTest && cY == yTest){
						cout << "call " << x << " " << y << endl;
					}
                    //printing = true;
					bool intersect = ClosestIntersection(cameraPos + (a-(lenseD-1)/2)* lenseSpacing * vec3( R[0][0], R[0][1], R[0][2]) + (b-(lenseD-1)/2)* lenseSpacing * vec3( R[1][0], R[1][1], R[1][2]) , glm::rotate(glm::rotate(
					glm::normalize(
					vec3( R[0][0], R[0][1], R[0][2] )*((float)(x-SCREEN_WIDTH/(2)))
					 + vec3( R[1][0], R[1][1], R[1][2] )*(float)((y-SCREEN_HEIGHT/(2))) + vec3( R[2][0], R[2][1], R[2][2] )*focalLength)
					, glm::atan((a-(lenseD-1)/2)*lenseSpacing/focalPoint),vec3(R[1][0], R[1][1], R[1][2])), glm::atan((a-(lenseD-1)/2)*lenseSpacing/focalPoint), vec3( R[0][0], R[0][1], R[0][2])), triangles, closestIntersection);

                    //printing = false;

                    vec3 color(0,0,0);
					if (intersect)
					{
						//printing = true;
						color = SoftLight(closestIntersection);
						//printing = false;
					}
					averageC += color;
				}
			}

			im[x][y][a*lenseD + b] = averageC/nsq;
			//PutPixelSDL( screen, x, y, averageC/nsq );
		}
	}

		}

	}

	for (int x = 0; x < SCREEN_WIDTH; x++)
	{
		for (int y = x % 2; y < SCREEN_HEIGHT; y+=2)
		{
			for (int z = 1;z<lenseD*lenseD;z++)
			{
				im[x][y][0] += im[x][y][z];
			}
			im[x][y][0] /= lenseD*lenseD;
			//cout << "(" << im[x][y][0].x << ", " << im[x][y][0].y << ", " << im[x][y][0].z << ")" << endl;
			//PutPixelSDL( screen, x, y, im[x][y][0] );
		}
	}

	for (int x = 0; x < SCREEN_WIDTH; x++)
	{
		for (int y = (x+1) % 2; y < SCREEN_HEIGHT; y+=2)
		{
			vec3 color = vec3(0,0,0);
			int count = 0;

			if (x>0)
			{
				color += im[x-1][y][0];
				count+=1;
			}
			if (y>0)
			{
				color += im[x][y-1][0];
				count+=1;
			}
			if (x<SCREEN_WIDTH-1)
			{
				color += im[x+1][y][0];
				count+=1;
			}
			if (y<SCREEN_HEIGHT-1)
			{
				color += im[x][y+1][0];
				count+=1;
			}
			im[x][y][0] = color/(float)count;
			//PutPixelSDL( screen, x, y, color/(float)count );
		}
	}

	for (int i = 0;i<SCREEN_WIDTH;i++){
		for (int j = 0;j<SCREEN_HEIGHT;j++){
			for (int h = 0;h<scale;h++){
				for (int v = 0;v<scale;v++){
					PutPixelSDL( screen, i*scale+h, j*scale+v, im[i][j][0]);
				}
			}
		}
	}


	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
