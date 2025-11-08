#ifndef GRAPHICS_H
#define GRAPHICS_H
#include <raylib.h>
#include "particles.h"

void initGraphics(int WIDTH, int HEIGHT);
void updateGraphics(Entities *simulationEntities, int screenDimensions[]);
#endif