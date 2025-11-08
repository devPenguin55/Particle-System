#include <stdio.h>
#include <stdlib.h>
#include <raylib.h>
#include <string.h>
#include "raymath.h"
#include "particles.h"

#define GRAVITY (Vector2){0, 300}
#define RADIUS 3

void initParticles(Entities *simulationEntities) {
    simulationEntities->capacity = 10;
    simulationEntities->amtEntities = 0;
    simulationEntities->entities = malloc(sizeof(Particle) * simulationEntities->capacity);

    simulationEntities->amtLinks = 0;
}

void spawnParticles(Entities *simulationEntities, int amount, int screenDimensions[]) {
    for (int i = 0; i < amount; i++) {
        if (simulationEntities->amtEntities >= simulationEntities->capacity) {
            simulationEntities->capacity *= 2;
            simulationEntities->entities = realloc(simulationEntities->entities, sizeof(Particle) * simulationEntities->capacity);
        }

        Particle *newParticle = &(simulationEntities->entities[(simulationEntities->amtEntities)++]);
        newParticle->positionCurrent = (Vector2){
            (int)(screenDimensions[0] / 2),
            (int)(screenDimensions[1] * 0.25)
        };
        newParticle->positionOld = (Vector2){
            // newParticle->positionCurrent.x + (simulationEntities->amtEntities % 8) - 4,
            newParticle->positionCurrent.x+0.025 ,
            newParticle->positionCurrent.y
        };
        
        newParticle->radius = RADIUS;
        newParticle->isStatic = 0;
        newParticle->color = (Color){
            (unsigned char)(127 * (sin(0.3*simulationEntities->amtEntities + 0) + 1)),
            (unsigned char)(127 * (sin(0.3*simulationEntities->amtEntities + 2) + 1)),
            (unsigned char)(127 * (sin(0.3*simulationEntities->amtEntities + 4) + 1)),
            255
        };

    }
}

void createIndependentLink(Entities *simulationEntities, int distance, int index1, int index2) {
    simulationEntities->links[(simulationEntities->amtLinks)++] = (Link){
        index1,
        index2,
        distance
    };
}

void createConnectedLinkChain(Entities *simulationEntities, Vector2 startPos, Vector2 endPos, double concentration) {
    int currentIndex = simulationEntities->amtEntities;
    // int initialIndex = simulationEntities->amtEntities;



    Vector2 direction = Vector2Subtract(startPos, endPos);
    float directionLength = Vector2Length(direction);

    Vector2 normal = Vector2Scale(direction, 1.0/directionLength);

    // direction (offset) then find the max amt of particles that could fit
    
    int amtParticles = ceil(concentration*(directionLength / (2*RADIUS+1)));
    float distanceBetweenParticles = 1+(directionLength - amtParticles*RADIUS*2) / ((amtParticles > 1) ? (amtParticles - 1) : 1);
    if (concentration == 1) {
        distanceBetweenParticles = 1;
    }
    for (int i=0; i<amtParticles; i++) {
        if (simulationEntities->amtEntities >= simulationEntities->capacity) {
            simulationEntities->capacity *= 2;
            simulationEntities->entities = realloc(simulationEntities->entities, sizeof(Particle) * simulationEntities->capacity);
        }
        Particle *newParticle = &(simulationEntities->entities[(simulationEntities->amtEntities)++]);
        newParticle->positionCurrent = Vector2Add(startPos, Vector2Scale(normal, -(i)*(RADIUS*2+distanceBetweenParticles)));

        newParticle->positionOld = (Vector2){
            newParticle->positionCurrent.x,
            newParticle->positionCurrent.y
        };
        newParticle->radius = RADIUS;
        newParticle->isStatic = 0;
        newParticle->color = (Color){
            (unsigned char)(127 * (sin(0.3*simulationEntities->amtEntities + 0) + 1)),
            (unsigned char)(127 * (sin(0.3*simulationEntities->amtEntities + 2) + 1)),
            (unsigned char)(127 * (sin(0.3*simulationEntities->amtEntities + 4) + 1)),
            255
        };
    }

    simulationEntities->entities[currentIndex].isStatic = 1;

    // 2 particles -> 1 chain
    // 3 particles -> 2 chains
    // 4 particles -> 3 chains
    // n particles -> n-1 chains
    for (int i = 0; i < (amtParticles-1); i++) {
        simulationEntities->links[(simulationEntities->amtLinks)++] = (Link){
            currentIndex,
            currentIndex+1,
            distanceBetweenParticles
        };
        currentIndex++;
    }

    // if (isHoriz) {
    //     simulationEntities->entities[simulationEntities->amtEntities-1].isStatic = 1;
    //     simulationEntities->entities[simulationEntities->amtEntities-1].positionCurrent.x += distance*(amtParticles-1)+(amtParticles-1)*2*(simulationEntities->entities[simulationEntities->amtEntities-1].radius);
    //     simulationEntities->entities[simulationEntities->amtEntities-1].positionOld.x += distance*(amtParticles-1)+(amtParticles-1)*2*(simulationEntities->entities[simulationEntities->amtEntities-1].radius);

    //     double halfDiff = 0.5*(simulationEntities->entities[simulationEntities->amtEntities-1].positionCurrent.x - simulationEntities->entities[initialIndex].positionCurrent.x);
    //     simulationEntities->entities[initialIndex].positionCurrent.x = screenDimensions[0] * 0.5 - halfDiff;
    //     simulationEntities->entities[initialIndex].positionOld.x = screenDimensions[0] * 0.5 - halfDiff;

    //     simulationEntities->entities[simulationEntities->amtEntities-1].positionCurrent.x = screenDimensions[0] * 0.5 + halfDiff;
    //     simulationEntities->entities[simulationEntities->amtEntities-1].positionOld.x = screenDimensions[0] * 0.5 + halfDiff;

    // }
}

void updateParticles(Entities *simulationEntities, double dt) {
    for (int particleIdx = 0; particleIdx < simulationEntities->amtEntities; particleIdx++) {
        Particle *curParticle = &(simulationEntities->entities[particleIdx]);
        if (curParticle->isStatic) {
            continue;
        }
        
        Vector2 velocity = Vector2Subtract(curParticle->positionCurrent, curParticle->positionOld);

        // velocity = Vector2Scale(velocity, 0.98);

        memcpy(&(curParticle->positionOld), &(curParticle->positionCurrent), sizeof(Vector2));

        curParticle->positionCurrent = Vector2Add(curParticle->positionCurrent, velocity);

        curParticle->positionCurrent = Vector2Add(
            curParticle->positionCurrent, 
            Vector2Scale(curParticle->acceleration, dt*dt)
        );

        
    }
}

void accelerateParticle(Particle *particle, Vector2 addAcceleration) {
    particle->acceleration = Vector2Add(particle->acceleration, addAcceleration);
}

void applyBoundaries(Entities *simulationEntities, int screenDimensions[], double dt) {
    Vector2 boundaryCenter = {screenDimensions[0] / 2, screenDimensions[1] / 2};
    double boundaryRadius = (screenDimensions[0] / 2);

    for (int particleIdx = 0; particleIdx < simulationEntities->amtEntities; particleIdx++) {
        Particle *curParticle = &(simulationEntities->entities[particleIdx]);
        if (curParticle->isStatic) {
            continue;
        }

        Vector2 offset = Vector2Subtract(curParticle->positionCurrent, boundaryCenter);

        // dist of particle to center      
        double dist = Vector2Length(offset) + curParticle->radius;
        // if dist is within radius, fine, if not, take it out
        if (dist > boundaryRadius) {
            // normalized direction -> the offset / dist
            Vector2 direction = Vector2Scale(offset, 1.0/dist);
            curParticle->positionCurrent = Vector2Add(boundaryCenter, Vector2Scale(direction, boundaryRadius));
        }  
    }
}

void applyCollisions(Entities *simulationEntities) {
    for (int particle1Idx = 0; particle1Idx < simulationEntities->amtEntities; particle1Idx++) {
        Particle *particle1 = &(simulationEntities->entities[particle1Idx]);

        for (int particle2Idx = particle1Idx+1; particle2Idx < simulationEntities->amtEntities; particle2Idx++) {
            Particle *particle2 = &(simulationEntities->entities[particle2Idx]);

            Vector2 offset = Vector2Subtract(particle2->positionCurrent, particle1->positionCurrent);
            // compute dist between centers
            double dist = Vector2Length(offset);
            if (dist < (particle1->radius + particle2->radius)) {
                // direction
                double penetration = (particle1->radius + particle2->radius) - dist;
                
                Vector2 normal;
                if (dist != 0.0) {
                    normal = Vector2Scale(offset, 1.0/dist);
                } else {
                    printf("had to do it!\n");
                    float jitter = 0.01f;
                    normal = (Vector2){ 
                        ((float)rand()/RAND_MAX - 0.5f) * jitter,
                        ((float)rand()/RAND_MAX - 0.5f) * jitter 
                    };
                }

                if (!particle1->isStatic) {
                    particle1->positionCurrent = Vector2Add(particle1->positionCurrent, Vector2Scale(normal, (-1*penetration)/2.0));
                }
                if (!particle2->isStatic) {
                    particle2->positionCurrent = Vector2Add(particle2->positionCurrent, Vector2Scale(normal, (1*penetration)/2.0));
                }
            }
        }
    }
}

void applyLinks(Entities *simulationEntities, int screenDimensions[]) {
    for (int linkIdx = 0; linkIdx < simulationEntities->amtLinks; linkIdx++) {
        Link *curLink = &(simulationEntities->links[linkIdx]);

        if (curLink->particle1Index == -1 && curLink->particle2Index == -1) {
            continue;
        }

        Particle *particle1 = &(simulationEntities->entities[curLink->particle1Index]);
        Particle *particle2 = &(simulationEntities->entities[curLink->particle2Index]);
        
        Vector2 offset = Vector2Subtract(particle1->positionCurrent, particle2->positionCurrent);
        double currentParticleDist = Vector2Length(offset);


        double delta = curLink->linkDistance - currentParticleDist + particle1->radius + particle2->radius;
        Vector2 normal = Vector2Scale(offset, 1.0/currentParticleDist);

        if (!particle1->isStatic) {
            particle1->positionCurrent = Vector2Add(particle1->positionCurrent, Vector2Scale(normal, 0.005*delta));
        }
        if (!particle2->isStatic) {
            particle2->positionCurrent = Vector2Add(particle2->positionCurrent, Vector2Scale(normal, -0.005*delta));
        }
    }
}

int circleIntersectsLine(Vector2 A, Vector2 B, Vector2 center, double radius) {
    Vector2 AB = {B.x - A.x, B.y - A.y};
    Vector2 AC = {center.x - A.x, center.y - A.y};

    float abLengthSquared = AB.x * AB.x + AB.y * AB.y;
    float t = (AC.x * AB.x + AC.y * AB.y) / abLengthSquared;

    if (t < 0) t = 0;
    else if (t > 1) t = 1;

    Vector2 closest = {A.x + AB.x * t, A.y + AB.y * t};

    float dx = closest.x - center.x;
    float dy = closest.y - center.y;
    float distSquared = dx * dx + dy * dy;

    return distSquared <= radius * radius;
}

void applyMouseForce(Entities *simulationEntities) {
    Vector2 mousePosition = GetMousePosition();
    for (int particleIdx = 0; particleIdx < simulationEntities->amtEntities; particleIdx++) {
        Particle *curParticle = &(simulationEntities->entities[particleIdx]);
        if (curParticle->isStatic) {
            continue;
        }

        Vector2 offset = Vector2Subtract(mousePosition, curParticle->positionCurrent);
        double dist = Vector2Length(offset);
        
        double mouseRadiusOfInfluence = 150;
        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) { 
            mouseRadiusOfInfluence = 25;
        }


        if (dist <= mouseRadiusOfInfluence) {
            // close enough to be considered
            double delta = mouseRadiusOfInfluence - dist;
            Vector2 normal = Vector2Scale(offset, 1.0/dist);

            if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
                // attract
                Vector2 vel = Vector2Subtract(curParticle->positionCurrent, curParticle->positionOld);
                vel = Vector2Add(Vector2Scale(vel, 0.9999), Vector2Scale(normal, 0.0002 * delta));
                curParticle->positionOld = Vector2Subtract(curParticle->positionCurrent, vel);
            } else if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
                // repel
                Vector2 vel = Vector2Subtract(curParticle->positionCurrent, curParticle->positionOld);
                vel = Vector2Add(vel, Vector2Scale(normal, -delta * 0.1)); 
                curParticle->positionOld = Vector2Subtract(curParticle->positionCurrent, vel);
            }
        }
    }

    if (IsMouseButtonDown(MOUSE_MIDDLE_BUTTON)) {
        for (int chainIdx = 0; chainIdx < simulationEntities->amtLinks; chainIdx++) {
            double mouseRadius = 5;

            if (circleIntersectsLine(
                simulationEntities->entities[simulationEntities->links[chainIdx].particle1Index].positionCurrent,
                simulationEntities->entities[simulationEntities->links[chainIdx].particle2Index].positionCurrent,
                mousePosition,
                mouseRadius)
            ) {
                // the mouse intersects with this link
                simulationEntities->links[chainIdx].particle1Index = -1;
                simulationEntities->links[chainIdx].particle2Index = -1;
            }
        }
    }
}

void solveState(Entities *simulationEntities, int screenDimensions[], double dt) {
    // return;
    int amtSubsteps = 8;
    double substepDt = dt / amtSubsteps;
    for (int substep = 0; substep < amtSubsteps; substep++) {
        for (int particleIdx = 0; particleIdx < simulationEntities->amtEntities; particleIdx++) {
            Particle *curParticle = &(simulationEntities->entities[particleIdx]);
            curParticle->acceleration.x = 0;
            curParticle->acceleration.y = 0;
    
            accelerateParticle(curParticle, GRAVITY);
        }
    
        updateParticles(simulationEntities, substepDt);
        applyBoundaries(simulationEntities, screenDimensions, substepDt);
        applyCollisions(simulationEntities);    
        applyLinks(simulationEntities, screenDimensions);
        applyMouseForce(simulationEntities);
    }
} 