/* 
  Object3d Example File 
  Displays a Spinning Icosahedron
*/

#include "SmartMatrix.h"
#include "FastLED.h"
#include "Object3d.h"
#include "icosahedron.h" 

#define kMatrixWidth  32
#define kMatrixHeight 32


#define NUM_LEDS (kMatrixWidth * kMatrixHeight)
#define DEFAULT_BRIGHTNESS 255

CRGB leds[kMatrixWidth * kMatrixHeight];

#define MESH_VERTICES 12
#define MESH_FACES 20
#define MESH_FOV 45
#define XCENTER 15.5
#define YCENTER 15.5

Object3d Icosa(MESH_VERTICES, MESH_FACES);

bool MESH_WIREFRAME = false;

float AngX = 20.0;
float AngY = 10.0;
uint32_t nextFrame = 0;
int cubeHue = 0;

void transformMesh() { 
  Icosa.Rotate(AngX,AngY,0);
  Icosa.Translate(0,0,-30);
  Icosa.Scale(1);
  Icosa.ApplyTransforms();
  Icosa.sortDepthMap();
  AngX +=.05;
  AngY +=.05;
}

void drawMesh(int Hue) { 
  int ax, ay,bx,by,cy,cx,V,baseHue;
  baseHue = Hue;
  for(int i=0; i<MESH_FACES; i++) {
    ax = Icosa.Render(XCENTER,MESH_FOV,i,"ax");
    ay = Icosa.Render(YCENTER,MESH_FOV,i,"ay");

    bx = Icosa.Render(XCENTER,MESH_FOV,i,"bx");
    by = Icosa.Render(YCENTER,MESH_FOV,i,"by");

    cx = Icosa.Render(XCENTER,MESH_FOV,i,"cx");
    cy = Icosa.Render(YCENTER,MESH_FOV,i,"cy");

    if(MESH_WIREFRAME) { 
      V = map(Icosa.depthMap[i].depth,Icosa.minDepth,Icosa.maxDepth,100,255);
      // replace this next line with your triangle drawing function
        pSmartMatrix->drawTriangle(ax,ay,bx,by,cx,cy, CRGB(CHSV(baseHue,255,V)));
        baseHue +=10;
    } else { 
      V = map(Icosa.depthMap[i].depth,Icosa.minDepth,Icosa.maxDepth,-100,255); // larger "range" as back faces won't be drawn
      V = constrain(V,0,255); // prevents flashing edges.
      // replace this next line with your triangle drawing function 
      pSmartMatrix->fillTriangle(ax,ay,bx,by,cx,cy, CRGB(CHSV(baseHue,255,V)));
      baseHue +=10;
    }      
  }
}

void drawIcosahedron() { 
  // clear the screen
  pSmartMatrix->fillScreen(CRGB(0,0,0));
  transformMesh();
  drawMesh(cubeHue);
  cubeHue++;
  nextFrame = millis() + 10;
}

void setup() { 
  LEDS.addLeds<SmartMatrix>(leds,NUM_LEDS);
  LEDS.setBrightness(DEFAULT_BRIGHTNESS);
  pSmartMatrix->setColorCorrection(cc24);
  Icosa.loadMesh(icosahedron_v, icosahedron_f); 
}

void loop() { 
  if(millis() >= nextFrame) { 
    drawIcosahedron();
    LEDS.show();
  }
}