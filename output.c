#include "output.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "watchtowerStruct.c"
#include "dcel.h"

#define NOFACE (-1)

void outputResult(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel){
    FILE *outputfile = fopen(outputfileName, "w");
    assert(outputfile);
    int i, j;
    int population = 0;
    if(!dcel){
        /* Simple path - avoids DCEL entirely. */
        for(i = 0; i < wtCount; i++){
            if(! wts[i]){
                continue;
            }
            fprintf(outputfile,
                "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
                "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
                (wts[i])->watchtowerID, 
                (wts[i])->postcode, 
                (wts[i])->populationServed, 
                (wts[i])->contact, 
                (wts[i])->x, 
                (wts[i])->y);
            population += wts[i]->populationServed;
        }
        fprintf(outputfile,
            "Face undefined population served: %d\n", population);
    } else {
        int faceCount = getFaceCount(dcel);
        int *faceMembership = (int *) malloc(sizeof(int)*wtCount);
        assert(faceMembership);
        for(i = 0; i < wtCount; i++){
            faceMembership[i] = NOFACE;
            for(j = 0; j < faceCount; j++){
                if(inFace(dcel, wts[i]->x, wts[i]->y, j)){
                    // if(faceMembership[i] != NOFACE){
                    //     printf("WT %d on edge of %d and %d!\n", i, faceMembership[i], j);
                    // }
                    faceMembership[i] = j;
                    break;
                }
            }
        }
        int *faceSize = (int *) malloc(sizeof(int)*(faceCount + 1));
        assert(faceSize);
        for(i = 0; i < (faceCount + 1); i++){
            faceSize[i] = 0;
        }
        for(i = 0; i < wtCount; i++){
            if(faceMembership[i] == NOFACE){
                faceSize[faceCount] = faceSize[faceCount] + 1;
            } else {
                faceSize[faceMembership[i]] = faceSize[faceMembership[i]] + 1;
            }
        }
        int **faceData = (int **) malloc(sizeof(int *)*(faceCount + 1));
        assert(faceData);
        for(i = 0; i < (faceCount + 1); i++){
            faceData[i] = (int *) malloc(sizeof(int)*faceSize[i]);
            assert(faceData[i] || faceSize[i] == 0);
            /* Reset the size of arrays so we can put them back in. */
            faceSize[i] = 0;
        }
        /* Put them into arrays. */
        for(i = 0; i < wtCount; i++){
            if(faceMembership[i] == NOFACE){
                faceData[faceCount][faceSize[faceCount]] = i;
                faceSize[faceCount] = faceSize[faceCount] + 1;
            } else {
                faceData[faceMembership[i]][faceSize[faceMembership[i]]] = i;
                faceSize[faceMembership[i]] = faceSize[faceMembership[i]] + 1;
            }
        }
        if(faceMembership){
            free(faceMembership);
        }
        /* Output */
        int *populationTotals = (int *) malloc(sizeof(int) * (faceCount + 1));
        assert(populationTotals);
        for(i = 0; i < (faceCount + 1); i++){
            populationTotals[i] = 0;
            // if(i != faceCount || faceSize[i] != 0){
            if(i != faceCount){
                fprintf(outputfile, "%d\n", i);
                for(j = 0; j < faceSize[i]; j++){
                    fprintf(outputfile,
                        "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
                        "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
                        (wts[faceData[i][j]])->watchtowerID, 
                        (wts[faceData[i][j]])->postcode, 
                        (wts[faceData[i][j]])->populationServed, 
                        (wts[faceData[i][j]])->contact, 
                        (wts[faceData[i][j]])->x, 
                        (wts[faceData[i][j]])->y);
                    populationTotals[i] += wts[faceData[i][j]]->populationServed;
                }
            }
        }
        
        for(i = 0; i < (faceCount); i++){
            fprintf(outputfile,
                "Face %d population served: %d\n", i, populationTotals[i]);
        }
        // if(faceSize[faceCount] > 0){
        //     fprintf(outputfile,
        //         "Face undefined population served: %d\n", populationTotals[faceCount]);
        // }
        if(populationTotals){
            free(populationTotals);
        }
        if(faceData){
            for(i = 0; i < (faceCount + 1); i++){
                if(faceData[i]){
                    free(faceData[i]);
                }
            }
            free(faceData);
        }
        if(faceSize){
            free(faceSize);
        }
    }
    fclose(outputfile);
}

void outputResultDiameter(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel){
    FILE *outputfile = fopen(outputfileName, "w");
    assert(outputfile);
    int i;
    
    /* Must have DCEL. */
    assert(dcel);
    
    for(i = 0; i < wtCount; i++){
        if(! wts[i]){
            continue;
        }
        // get diameter and store in watchtower structs
        wts[i]->diameter = getDiameter(dcel, wts[i]->face);
        fprintf(outputfile,
            "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
            "Watchtower Point of Contact Name: %s, x: %f, y: %f, Diameter of Cell: %f\n",
            (wts[i])->watchtowerID, 
            (wts[i])->postcode, 
            (wts[i])->populationServed, 
            (wts[i])->contact, 
            (wts[i])->x, 
            (wts[i])->y,
            (wts[i])->diameter);
    }
    
    fclose(outputfile);
}

void outputResultDiameterSorted(char *outputfileName, struct watchtowerStruct **wts, int wtCount, 
    struct DCEL *dcel){

    int i, j;

    FILE *outputfile = fopen(outputfileName, "w");
    assert(outputfile);
    
    /* Must have DCEL. */
    assert(dcel);

    // get diameter and store in watchtower structs
    for(i = 0; i < wtCount; i++){
        if(! wts[i]){
            continue;
        }
        wts[i]->diameter = getDiameter(dcel, wts[i]->face);
    }

    // insertion sort by diamater
    for (i=1; i<wtCount; i++) {
        for (j=i-1; j>=0 && wts[j+1]->diameter < wts[j]->diameter; j--) {
            wtsSwap(&wts[j+1], &wts[j]);
        }
    }

    // print watchtower data to file in sorted order
    for(i = 0; i < wtCount; i++){
        if(! wts[i]){
            continue;
        }
        fprintf(outputfile,
                "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
                "Watchtower Point of Contact Name: %s, x: %f, y: %f, Diameter of Cell: %f\n",
                (wts[i])->watchtowerID,
                (wts[i])->postcode,
                (wts[i])->populationServed,
                (wts[i])->contact,
                (wts[i])->x,
                (wts[i])->y,
                (wts[i])->diameter);
    }
    
    fclose(outputfile);
}

void wtsSwap(struct watchtowerStruct **wts1, struct watchtowerStruct **wts2) {
    struct watchtowerStruct *temp = *wts1;
    *wts1 = *wts2;
    *wts2 = temp;
}

char *getWTDataString(struct watchtowerStruct **wts, int wtIndex){
    char *s = NULL;
    /* Get size */
    size_t requiredSpace = snprintf(NULL, 0, 
        "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
        "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
        (wts[wtIndex])->watchtowerID, 
        (wts[wtIndex])->postcode, 
        (wts[wtIndex])->populationServed, 
        (wts[wtIndex])->contact, 
        (wts[wtIndex])->x, 
        (wts[wtIndex])->y);
    s = (char *) malloc(sizeof(char) * (requiredSpace + 1));
    assert(s);
    sprintf(s, 
        "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
        "Watchtower Point of Contact Name: %s, x: %f, y: %f\n",
        (wts[wtIndex])->watchtowerID, 
        (wts[wtIndex])->postcode, 
        (wts[wtIndex])->populationServed, 
        (wts[wtIndex])->contact, 
        (wts[wtIndex])->x, 
        (wts[wtIndex])->y);
    
    return s;
}

void freeWTDataString(char *s){
    if(s){
        free(s);
    }
}

double getWTx(int index, struct watchtowerStruct **wts){
    if(!wts || !(wts[index])){
        return 0;
    }
    
    return (wts[index])->x;
}

double getWTy(int index, struct watchtowerStruct **wts){
    if(!wts || !(wts[index])){
        return 0;
    }
    
    return (wts[index])->y;
}
