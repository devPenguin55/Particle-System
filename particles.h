#ifndef PARTICLES_H
#define PARTICLES_H

typedef struct Particle
{
    Vector2 positionCurrent;
    Vector2 positionOld;
    Vector2 acceleration;
    int radius;
    int isStatic;
    Color color;
} Particle;

typedef struct Link
{
    int particle1Index;
    int particle2Index;
    int linkDistance;
} Link;


typedef struct Entities {
    Particle *entities;
    Link links[100000];
    int amtEntities;
    int capacity;
    int amtLinks;
} Entities;

void initParticles(Entities *simulationEntities);
void spawnParticles(Entities *simulationEntities, int amount, int screenDimensions[]);
void createIndependentLink(Entities *simulationEntities, int distance, int index1, int index2);
void createConnectedLinkChain(Entities *simulationEntities, Vector2 startPos, Vector2 endPos, double concentration);
void updateParticles(Entities *simulationEntities, double dt);
int circleIntersectsLine(Vector2 A, Vector2 B, Vector2 center, double radius);
void applyMouseForce(Entities *simulationEntities);
void solveState(Entities *simulationEntities, int screenDimensions[], double dt);
void applyCollisions(Entities *simulationEntities);

#endif