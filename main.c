#include <stdio.h>
#include <stdlib.h>
#include "graphics.h"
#include "particles.h"
#include <math.h>

int main() {
    int screenDimensions[] = {800, 800};
    initGraphics(screenDimensions[0], screenDimensions[1]);

    Entities simulationEntities;
    initParticles(&simulationEntities);

    // spawnParticles(&simulationEntities, 4, screenDimensions);

    // int startIndex = 0;
    // float x0 = 100+startIndex*150, y0 = 100; 
    // float s = 30;       
    
    // createIndependentLink(&simulationEntities, 35, 0+startIndex, 1+startIndex); 
    // createIndependentLink(&simulationEntities, 35, 1+startIndex, 3+startIndex); 
    // createIndependentLink(&simulationEntities, 35, 3+startIndex, 2+startIndex);
    // createIndependentLink(&simulationEntities, 35, 2+startIndex, 0+startIndex);
    // createIndependentLink(&simulationEntities, 35, 0+startIndex, 3+startIndex); 
    // createIndependentLink(&simulationEntities, 35, 1+startIndex, 2+startIndex); 
    // simulationEntities.entities[0].positionCurrent = (Vector2){x0, y0};
    // simulationEntities.entities[1].positionCurrent = (Vector2){x0 + s, y0};
    // simulationEntities.entities[2].positionCurrent = (Vector2){x0, y0 + s};
    // simulationEntities.entities[3].positionCurrent = (Vector2){x0 + s, y0 + s};
    // simulationEntities.entities[0].positionOld = (Vector2){x0, y0};
    // simulationEntities.entities[1].positionOld = (Vector2){x0 + s, y0};
    // simulationEntities.entities[2].positionOld = (Vector2){x0, y0 + s};
    // simulationEntities.entities[3].positionOld = (Vector2){x0 + s, y0 + s};

    // createConnectedLinkChain(&simulationEntities, (Vector2){350, 40}, (Vector2){350, 250}, 1);
    
    
    
    
    spawnParticles(&simulationEntities, 5000, screenDimensions);
    
    // int ROWS = 30;
    // int COLS = 17;
    // spawnParticles(&simulationEntities, ROWS*COLS, screenDimensions);
    
    // float x0 = 160; 
    // float y0 = 100; 
    // float s = 20; 
    // for (int row = 0; row < ROWS; row++) {
    //     for (int col = 0; col < COLS; col++) {
    //         int index = row * COLS + col;

    //         printf("x %d y %d index %d\n", col, row, index);

    //         if (index < COLS) {
    //             simulationEntities.entities[index].isStatic = 1;
    //         }

    //         simulationEntities.entities[index].positionCurrent = (Vector2){x0 + col * (s+10), y0 + row * (s+10)};
    //         simulationEntities.entities[index].positionOld     = (Vector2){x0 + col * (s+10), y0 + row * (s+10)};

    //         if (col < (COLS-1)) {
    //             createIndependentLink(&simulationEntities, s, index, index + 1);
    //         }

    //         if (row < (ROWS-1)) {
    //             createIndependentLink(&simulationEntities, s, index, index + COLS);
    //         }

    //         if ((col < (COLS-1)) && (row < (ROWS-1))) {
    //             createIndependentLink(&simulationEntities, sqrtf(s*s + s*s), index, index + COLS + 1);
    //         }

    //         if ((col > 0) && (row < (ROWS-1))) {
                
    //             createIndependentLink(&simulationEntities, sqrtf(s*s + s*s), index, index + COLS - 1);
    //         }
    //     }
    // }



    int tick = 0;
    while (!WindowShouldClose()) {
        double dt = GetFrameTime();
        solveState(&simulationEntities, screenDimensions, dt);
        updateGraphics(&simulationEntities, screenDimensions);
        
        // if ((tick % 10) == 0) {
        //     spawnParticles(&simulationEntities, 1, screenDimensions);
        // }
        tick++;

    }

    CloseWindow();

    return 0;
}