#include <stdio.h>
#include <stdlib.h>
#include "graphics.h"
#include "particles.h"

void initGraphics(int WIDTH, int HEIGHT) {
    InitWindow(WIDTH, HEIGHT+20, "Particle Simulation");
    SetTargetFPS(60);
}

void updateGraphics(Entities *simulationEntities, int screenDimensions[]) {
    BeginDrawing();
    ClearBackground(BLACK);

    DrawCircleV((Vector2){screenDimensions[0] / 2, screenDimensions[1] / 2}, (screenDimensions[0] / 2), DARKGRAY);

    for (int particleIdx = 0; particleIdx < simulationEntities->amtEntities; particleIdx++) {
        Particle *curParticle = &(simulationEntities->entities[particleIdx]);
        DrawCircleV(curParticle->positionCurrent, curParticle->radius, curParticle->color);
        // DrawCircleV(curParticle->positionCurrent, curParticle->radius, LIGHTGRAY);
    }

    for (int chainIdx = 0; chainIdx < simulationEntities->amtLinks; chainIdx++) {
        DrawLineV(
            simulationEntities->entities[simulationEntities->links[chainIdx].particle1Index].positionCurrent,
            simulationEntities->entities[simulationEntities->links[chainIdx].particle2Index].positionCurrent,
            WHITE
        );
    }

    
    EndDrawing();
}