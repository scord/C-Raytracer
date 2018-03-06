#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>
#include <SDL.h>

// Used to describe a triangular surface:
class Triangle
{
public:
	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 v2;
	glm::vec3 normal;
	glm::vec3 color;
	glm::vec3** texture;
	int texWidth;
	int texHeight;
	float spec;
	float diff;
	float shin;
	float refl;
	bool textured;


	Triangle( glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 color, float diff, float spec, float shin, float refl, int texWidth, int texHeight, glm::vec3** texture)
		: v0(v0), v1(v1), v2(v2), color(color), diff(diff), spec(spec), refl(refl), shin(shin), texWidth(texWidth), texHeight(texHeight), texture(texture)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
		glm::vec3 e1 = v1-v0;
		glm::vec3 e2 = v2-v0;
		normal = glm::normalize( glm::cross( e2, e1) );
	}
};

glm::vec3** getTexture(const char* filename)
{
	using glm::vec3;

	SDL_Surface* surface;
	surface = SDL_LoadBMP(filename);
	Uint32 *pixels = (Uint32 *)surface->pixels;

	vec3** texture;
	texture = new vec3 * [surface->w];



	for (int i = 0; i < surface->w; i++)
	{
		texture[i] = new vec3 [surface->h];
		assert(texture[i]);
	}

	for (int x = 0; x < surface->w; x++)
	{
		for (int y = 0; y < surface->h; y++)
		{
			SDL_Color col;

			Uint32 pixel = pixels[(y*surface->w)+x];

			SDL_GetRGB(pixel, surface->format, &col.r, &col.g, &col.b);
			Uint8 *colors = (Uint8*)&pixel;
			texture[x][y] = vec3((float)col.r/255.0f, (float)col.g/255.0f, (float)col.b/255.0f);
		}
	}
	return texture;
}

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles )
{
	using glm::vec3;

	SDL_Init(SDL_INIT_VIDEO);

	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );
	vec3 black( 0.0f, 0.0f, 0.0f);

	triangles.clear();
	triangles.reserve( 5*2*3 );

	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	vec3 A(L,0,0);
	vec3 B(0,0,0);
	vec3 C(L,0,L);
	vec3 D(0,0,L);

	vec3 E(L,L,0);
	vec3 F(0,L,0);
	vec3 G(L,L,L);
	vec3 H(0,L,L);

	vec3** texture = getTexture("concrete.bmp");
	vec3** brickTexture = getTexture("bricks.bmp");

	std::cout << sizeof(texture) << std::endl;
	// Floor:
	triangles.push_back( Triangle( C, B, A, white, 0.1, 0.9, 100, 0, 800, 800, texture) );
	triangles.push_back( Triangle( C, D, B, white, 0.1, 0.9, 100, 0, 800, 800, texture) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, white, 1, 0, 0, 0, 800, 800, texture) );
	triangles.push_back( Triangle( C, E, G, white, 1, 0, 0, 0, 800, 800, texture) );

	// Right wall
	triangles.push_back( Triangle( F, B, D, white, 1, 0, 0, 0, 800, 800, texture) );
	triangles.push_back( Triangle( H, F, D, white, 1, 0, 0, 0, 800, 800, texture) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G, white, 1, 0, 0, 0, 800, 800, texture) );
	triangles.push_back( Triangle( F, H, G, white, 1, 0, 0, 0, 800, 800, texture) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, white, 0, 1, 10, 0, 800, 800, texture) );
	triangles.push_back( Triangle( G, H, D, white, 0, 1, 10, 0, 800, 800, texture) );

	// ---------------------------------------------------------------------------
	// Short block

	A = vec3(290,0,114);
	B = vec3(130,0, 65);
	C = vec3(240,0,272);
	D = vec3( 82,0,225);

	E = vec3(290,165,114);
	F = vec3(130,165, 65);
	G = vec3(240,165,272);
	H = vec3( 82,165,225);

	// Front
	triangles.push_back( Triangle(E,B,A,red, 0.5, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(E,F,B,red, 0.5, 0.5, 5, 0, 800, 800, texture) );

	// Front
	triangles.push_back( Triangle(F,D,B,red, 0.5, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(F,H,D,red, 0.5, 0.5, 5, 0, 800, 800, texture) );

	// BACK
	triangles.push_back( Triangle(H,C,D,red, 0.5, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(H,G,C,red, 0.5, 0.5, 5, 0, 800, 800, texture) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,red, 0.5, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(E,A,C,red, 0.5, 0.5, 5, 0, 800, 800, texture) );

	// TOP
	triangles.push_back( Triangle(G,F,E,red, 0.5, 0.5, 5, 1, 800, 800, texture) );
	triangles.push_back( Triangle(G,H,F,red, 0.5, 0.5, 5, 1, 800, 800, texture) );

	// ---------------------------------------------------------------------------
	// Tall block

	A = vec3(423,0,247);
	B = vec3(265,0,296);
	C = vec3(472,0,406);
	D = vec3(314,0,456);

	E = vec3(423,330,247);
	F = vec3(265,330,296);
	G = vec3(472,330,406);
	H = vec3(314,330,456);

	// Front
	triangles.push_back( Triangle(E,B,A,blue, 0.1, 0.9, 5, 0.9, 800, 800, texture) );
	triangles.push_back( Triangle(E,F,B,blue, 0.1, 0.9, 5, 0.9, 800, 800, texture) );

	// Front
	triangles.push_back( Triangle(F,D,B,blue, 1, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(F,H,D,blue, 1, 0.5, 5, 0, 800, 800, texture) );

	// BACK
	triangles.push_back( Triangle(H,C,D,blue, 1, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(H,G,C,blue, 1, 0.5, 5, 0, 800, 800, texture) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,blue, 1, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(E,A,C,blue, 1, 0.5, 5, 0, 800, 800, texture) );

	// TOP
	triangles.push_back( Triangle(G,F,E,blue, 1, 0.5, 5, 0, 800, 800, texture) );
	triangles.push_back( Triangle(G,H,F,blue, 1, 0.5, 5, 0, 800, 800, texture) );


	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec3(1,1,1);
		triangles[i].v1 -= vec3(1,1,1);
		triangles[i].v2 -= vec3(1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].ComputeNormal();
	}
}

#endif
