#include "cubiomes/finders.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>
#include <string.h> // For memcpy

// -----------------------------------------------------------------------------
// Global configuration / defaults
#define MAX_SEEDS_TO_FIND 1 
int seedsFound = 0;          // Tracks how many seeds have been found so far

// Range for seed scanning:
uint64_t starting_seed = 1;
uint64_t end_seed      = 0ULL;  // Will be set in main (starting_seed + some offset)

// Search radius for scanning
int searchRadius = 1000;  // e.g., +/- 1000 around spawn
int useSpawn     = 1;     // 1 = use real spawn, 0 = use (customX, customZ)
int customX      = 0;
int customZ      = 0;

// Number of threads (can be set by user)
int tasksCount = 1;

// -----------------------------------------------------------------------------
// Union–find utility for clustering
int findSet(int parent[], int i)
{
    if (parent[i] != i)
        parent[i] = findSet(parent, parent[i]);
    return parent[i];
}
void unionSets(int parent[], int x, int y)
{
    int rx = findSet(parent, x);
    int ry = findSet(parent, y);
    if (rx != ry) parent[ry] = rx;
}

// -----------------------------------------------------------------------------
// Biome ID -> Name
const char* getBiomeName(int id)
{
    switch(id) {
        case 0:   return "Ocean";
        case 1:   return "Plains";
        case 2:   return "Desert";
        case 3:   return "Windswept Hills";
        case 4:   return "Forest";
        case 5:   return "Taiga";
        case 6:   return "Swamp";
        case 7:   return "River";
        case 10:  return "Frozen Ocean";
        case 11:  return "Frozen River";
        case 12:  return "Snowy Plains";
        case 13:  return "Snowy Mountains";
        case 14:  return "Mushroom Fields";
        case 15:  return "Mushroom Fields Shore";
        case 16:  return "Beach";
        case 17:  return "Desert Hills";
        case 18:  return "Windswept Forest";
        case 19:  return "Taiga Hills";
        case 20:  return "Mountain Edge";
        case 21:  return "Jungle";
        case 22:  return "Jungle Hills";
        case 23:  return "Sparse Jungle";
        case 24:  return "Deep Ocean";
        case 25:  return "Stony Shore";
        case 26:  return "Snowy Beach";
        case 27:  return "Birch Forest";
        case 28:  return "Birch Forest Hills";
        case 29:  return "Dark Forest";
        case 30:  return "Snowy Taiga";
        case 31:  return "Snowy Taiga Hills";
        case 32:  return "Old Growth Pine Taiga";
        case 33:  return "Giant Tree Taiga Hills";
        case 34:  return "Wooded Mountains";
        case 35:  return "Savanna";
        case 36:  return "Savanna Plateau";
        case 37:  return "Badlands";
        case 38:  return "Wooded Badlands";
        case 39:  return "Badlands Plateau";
        case 44:  return "Warm Ocean";
        case 45:  return "Lukewarm Ocean";
        case 46:  return "Cold Ocean";
        case 47:  return "Deep Warm Ocean";
        case 48:  return "Deep Lukewarm Ocean";
        case 49:  return "Deep Cold Ocean";
        case 50:  return "Deep Frozen Ocean";
        case 129: return "Sunflower Plains";
        case 130: return "Desert Lakes";
        case 131: return "Windswept Gravelly Hills";
        case 132: return "Flower Forest";
        case 133: return "Taiga Mountains";
        case 134: return "Swamp Hills";
        case 140: return "Ice Spikes";
        case 149: return "Modified Jungle";
        case 151: return "Modified Jungle Edge";
        case 155: return "Old Growth Birch Forest";
        case 156: return "Tall Birch Hills";
        case 157: return "Dark Forest Hills";
        case 158: return "Snowy Taiga Mountains";
        case 160: return "Old Growth Spruce Taiga";
        case 161: return "Giant Spruce Taiga Hills";
        case 162: return "Gravelly Mountains+";
        case 163: return "Windswept Savanna";
        case 164: return "Shattered Savanna Plateau";
        case 165: return "Eroded Badlands";
        case 166: return "Modified Wooded Badlands Plateau";
        case 167: return "Modified Badlands Plateau";
        case 168: return "Bamboo Jungle";
        case 169: return "Bamboo Jungle Hills";
        case 170: return "Underground";
        case 171: return "Underground";
        case 172: return "Underground";
        case 174: return "Dripstone Caves";
        case 175: return "Lush Caves";
        case 177: return "Meadow";
        case 178: return "Grove";
        case 179: return "Snowy Slopes";
        case 180: return "Frozen Peaks";
        case 181: return "Jagged Peaks";
        case 182: return "Stony Peaks";
        case 183: return "Deep Dark";
        case 184: return "Mangrove Swamp";
        case 185: return "Cherry Grove";
        case 186: return "Pale Garden";
        default:  return "Unknown Biome";
    }
}

// -----------------------------------------------------------------------------
// Approx. function to measure the "patch size" of a given biome around (x,z).
int getBiomePatchSize(Generator *g, int x, int z, int biome_id)
{
    int radius = 128; 
    int sx = (x - radius) >> 2;
    int sz = (z - radius) >> 2;
    int w = radius >> 1;
    int h = radius >> 1;
    int minX = INT_MAX, maxX = INT_MIN;
    int minZ = INT_MAX, maxZ = INT_MIN;
    int found = 0;

    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            int bx = sx + i;
            int bz = sz + j;
            int id = getBiomeAt(g, 4, bx, 0, bz);
            if (id == biome_id) {
                found = 1;
                if (bx < minX) minX = bx;
                if (bx > maxX) maxX = bx;
                if (bz < minZ) minZ = bz;
                if (bz > maxZ) maxZ = bz;
            }
            else if (found && (bx > maxX + 2 || bz > maxZ + 2)) {
                break;
            }
        }
    }
    if (!found) return 0;

    int sizeX = (maxX - minX + 1) * 4;
    int sizeZ = (maxZ - minZ + 1) * 4;
    return (sizeX + sizeZ) / 2;
}

// -----------------------------------------------------------------------------
// Required structure conditions
typedef struct {
    int structureType;   
    int minCount;
    int minHeight;
    int maxHeight;
    int requiredBiome;   
    int minBiomeSize;    
    int maxBiomeSize;    
} StructureRequirement;

// Per-biome size config for required patches
typedef struct {
    int biomeId;
    int minSize; // -1 if no minimum
    int maxSize; // -1 if no maximum
} BiomeSizeConfig;

// For "required" biomes
typedef struct {
    int *biomeIds;
    int biomeCount;
    BiomeSizeConfig *sizeConfigs;
    int configCount;
    int logCenters;
} BiomeRequirement;

// For "clustered" biomes
typedef struct {
    int *biomeIds;
    int biomeCount;
    int minSize;
    int maxSize;
    int logCenters;
} BiomeCluster;

// Top-level container for biome searches
typedef struct {
    BiomeRequirement *required;
    int requiredCount;
    BiomeCluster *clusters;
    int clusterCount;
} BiomeSearch;

// -----------------------------------------------------------------------------
// Example data for required/clustered biomes
static int reqGroup0[] = {};//{185}; // placeholder
static BiomeSizeConfig reqSizeConfigs[] = {{0,0,0}};//{{185, -1, -1}};
static BiomeRequirement reqGroup;

static int clusterGroup0[] = {};//{1, 185}; // e.g. {129, 185}
static const BiomeCluster clustGroup0 = {
    .biomeIds   = NULL, //clusterGroup0,
    .biomeCount = 0,//sizeof(clusterGroup0) / sizeof(clusterGroup0[0]),
    .minSize    = -1,//2,
    .maxSize    = -1,
    .logCenters = 1
};

// We'll store them as arrays
static BiomeRequirement requiredBiomes[1];
static int requiredBiomesCount = 1;

static const BiomeCluster biomeClusters[] = { clustGroup0 };
static const int biomeClustersCount = sizeof(biomeClusters)/sizeof(biomeClusters[0]);

static BiomeSearch biomeSearch = {
    .required      = requiredBiomes,
    .requiredCount = 1,
    .clusters      = (BiomeCluster *) biomeClusters,
    .clusterCount  = 1
};

// -----------------------------------------------------------------------------
// Structure cluster requirement
typedef struct {
    bool enabled;
    int clusterDistance;
    int *structureTypes;
    int count;
    int minClusterSize;  
} ClusterRequirement;

// Example: cluster these structure types
int clusterTypesArray[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 16, 18 };

ClusterRequirement clusterReq = {
    .enabled         = true,
    .clusterDistance = 32,
    .structureTypes  = clusterTypesArray,
    .count           = sizeof(clusterTypesArray)/sizeof(clusterTypesArray[0]),
    .minClusterSize  = 3
};

// -----------------------------------------------------------------------------
// Found positions
typedef struct {
    int structureType; 
    int x;
    int z;
} StructurePos;

// -----------------------------------------------------------------------------
// Optional invalid combination
typedef struct {
    int *types;
    int count;
} InvalidCombination;

#define NUM_INVALID_COMBINATIONS 2
int invComb1Arr[] = {16, 11}; 
int invComb2Arr[] = {16, 6};  
InvalidCombination invalidCombinations[NUM_INVALID_COMBINATIONS] = {
    { invComb1Arr, sizeof(invComb1Arr)/sizeof(invComb1Arr[0]) },
    { invComb2Arr, sizeof(invComb2Arr)/sizeof(invComb2Arr[0]) }
};

// Example list of structure requirements:
StructureRequirement structureRequirements[] = {
    //{ 5, 1, -10000, 10000, 1, -1, -1 }
};
int NUM_STRUCTURE_REQUIREMENTS = sizeof(structureRequirements)/sizeof(structureRequirements[0]);

// -----------------------------------------------------------------------------
// Helper compare function
int compareInts(const void *a, const void *b)
{
    int A = *(const int*)a;
    int B = *(const int*)b;
    return (A - B);
}

// -----------------------------------------------------------------------------
// Check if an invalid combination is a subset of the given cluster
bool isInvalidClusterDynamic(int *groupTypes, int groupSize)
{
    for (int ic = 0; ic < NUM_INVALID_COMBINATIONS; ic++) {
        int *invalidSet   = invalidCombinations[ic].types;
        int invalidCount  = invalidCombinations[ic].count;
        bool isSubset     = true;
        for (int k = 0; k < invalidCount; k++) {
            bool found = false;
            for (int j = 0; j < groupSize; j++) {
                if (groupTypes[j] == invalidSet[k]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                isSubset = false;
                break;
            }
        }
        if (isSubset)
            return true;
    }
    return false;
}

// -----------------------------------------------------------------------------
// scanBiomes: For required & clustered biome conditions (if any)
bool scanBiomes(Generator *g, int x0, int z0, int x1, int z1, BiomeSearch *bs)
{
    if (bs->requiredCount == 0 && bs->clusterCount == 0) {
        fprintf(stderr, "Error: No biome requirements provided.\n");
        return false;
    }
    bool success = true; // tracks if all required are found
    int step = 4;

    // 1) Required biome patches
    if (bs->requiredCount > 0) {
        bool anyRequiredFound = false;
        for (int i = 0; i < bs->requiredCount; i++) {
            BiomeRequirement *req = &bs->required[i];
            // Collect all cells that match any biome in req->biomeIds
            int capacity = 128, count = 0;
            StructurePos *positions = malloc(capacity * sizeof(StructurePos));
            if (!positions) { perror("malloc"); exit(1); }

            for (int zz = z0; zz <= z1; zz += step) {
                for (int xx = x0; xx <= x1; xx += step) {
                    int biome = getBiomeAt(g, 4, xx >> 2, 0, zz >> 2);
                    // check if biome is in req->biomeIds
                    for (int b = 0; b < req->biomeCount; b++) {
                        if (biome == req->biomeIds[b]) {
                            if (count == capacity) {
                                capacity *= 2;
                                positions = realloc(positions, capacity * sizeof(StructurePos));
                                if (!positions) { perror("realloc"); exit(1); }
                            }
                            positions[count].structureType = biome;
                            positions[count].x = xx;
                            positions[count].z = zz;
                            count++;
                            break; 
                        }
                    }
                }
            }

            if (count == 0) {
                // none found for this requirement
                free(positions);
                continue;
            }

            // union-find
            int *parent = malloc(count * sizeof(int));
            for (int c = 0; c < count; c++)
                parent[c] = c;

            for (int c = 0; c < count; c++) {
                for (int d = c + 1; d < count; d++) {
                    int dx = abs(positions[c].x - positions[d].x);
                    int dz = abs(positions[c].z - positions[d].z);
                    // same biome type, within "touching" distance => union
                    if (dx <= step && dz <= step &&
                        positions[c].structureType == positions[d].structureType)
                    {
                        unionSets(parent, c, d);
                    }
                }
            }

            bool *processed = calloc(count, sizeof(bool));
            for (int c = 0; c < count; c++) {
                int root = findSet(parent, c);
                if (processed[root]) continue;
                processed[root] = true;

                // gather stats
                double sumX = 0, sumZ = 0;
                int compCount = 0;
                int theBiome = positions[c].structureType;

                for (int e = 0; e < count; e++) {
                    if (findSet(parent, e) == root) {
                        sumX += positions[e].x;
                        sumZ += positions[e].z;
                        compCount++;
                    }
                }
                double centerX = sumX / compCount;
                double centerZ = sumZ / compCount;

                // check per-biome size constraints (if any)
                bool patchOk = true;
                // see if there's a config entry for this biome
                for (int sc = 0; sc < req->configCount; sc++) {
                    if (req->sizeConfigs[sc].biomeId == theBiome) {
                        int mn = req->sizeConfigs[sc].minSize;
                        int mx = req->sizeConfigs[sc].maxSize;
                        if (mn > -1 && compCount < mn) patchOk = false;
                        if (mx > -1 && compCount > mx) patchOk = false;
                        break;
                    }
                }

                if (patchOk) {
                    anyRequiredFound = true;
                    if (req->logCenters) {
                        printf("Required biome patch (reqIndex=%d, biome=%s): center (%.1f,%.1f), cells=%d\n",
                               i, getBiomeName(theBiome), centerX, centerZ, compCount);
                    }
                }
            }
            free(processed);
            free(parent);
            free(positions);
        } // end for each required group

        if (!anyRequiredFound) success = false;
    }

    // 2) Clustered biomes
    if (bs->clusterCount > 0) {
        for (int i = 0; i < bs->clusterCount; i++) {
            BiomeCluster *cl = &bs->clusters[i];
            // gather all cells that match any biome in cl->biomeIds
            int capacity = 128, count = 0;
            StructurePos *positions = malloc(capacity * sizeof(StructurePos));
            if (!positions) { perror("malloc"); exit(1); }

            for (int zz = z0; zz <= z1; zz += step) {
                for (int xx = x0; xx <= x1; xx += step) {
                    int biome = getBiomeAt(g, 4, xx >> 2, 0, zz >> 2);
                    for (int b = 0; b < cl->biomeCount; b++) {
                        if (biome == cl->biomeIds[b]) {
                            if (count == capacity) {
                                capacity *= 2;
                                positions = realloc(positions, capacity*sizeof(StructurePos));
                                if (!positions) { perror("realloc"); exit(1); }
                            }
                            positions[count].structureType = biome;
                            positions[count].x = xx;
                            positions[count].z = zz;
                            count++;
                            break; 
                        }
                    }
                }
            }

            // if none found, skip
            if (count == 0) {
                free(positions);
                continue;
            }

            // union-find them by adjacency
            int *parent = malloc(count * sizeof(int));
            for (int c = 0; c < count; c++)
                parent[c] = c;
            for (int c = 0; c < count; c++) {
                for (int d = c+1; d < count; d++) {
                    int dx = abs(positions[c].x - positions[d].x);
                    int dz = abs(positions[c].z - positions[d].z);
                    if (dx <= step && dz <= step) {
                        unionSets(parent, c, d);
                    }
                }
            }

            bool *processed = calloc(count, sizeof(bool));
            for (int c = 0; c < count; c++) {
                int root = findSet(parent, c);
                if (processed[root]) continue;
                processed[root] = true;

                double sumX = 0, sumZ = 0;
                int compCount = 0;
                bool usedBiome[256];
                memset(usedBiome, false, sizeof(usedBiome));

                for (int e = 0; e < count; e++) {
                    if (findSet(parent, e) == root) {
                        sumX += positions[e].x;
                        sumZ += positions[e].z;
                        compCount++;
                        int bId = positions[e].structureType;
                        if (bId >= 0 && bId < 256) {
                            usedBiome[bId] = true;
                        }
                    }
                }

                // how many distinct biomes?
                int distinctIDs[256];
                int dIdx = 0;
                for (int bId = 0; bId < 256; bId++) {
                    if (usedBiome[bId]) {
                        distinctIDs[dIdx++] = bId;
                    }
                }

                if (dIdx == 2) {
                    double cx = sumX / compCount;
                    double cz = sumZ / compCount;
                    bool sizeOk = true;
                    if (cl->minSize > -1 && compCount < cl->minSize) sizeOk = false;
                    if (cl->maxSize > -1 && compCount > cl->maxSize) sizeOk = false;

                    if (sizeOk && cl->logCenters) {
                        const char *bname1 = getBiomeName(distinctIDs[0]);
                        const char *bname2 = getBiomeName(distinctIDs[1]);
                        printf("Clustered biome group %d: %s + %s, center(%.1f,%.1f), count=%d\n",
                               i+1, bname1, bname2, cx, cz, compCount);
                    }
                }
            }
            free(processed);
            free(parent);
            free(positions);
        }
    }

    return success;
}

// -----------------------------------------------------------------------------
// Shared concurrency variables
volatile bool foundValidSeed    = false;
pthread_mutex_t seedMutex       = PTHREAD_MUTEX_INITIALIZER;
uint64_t validSeed              = 0;
volatile uint64_t currentSeed   = 0; // threads will increment this

// -----------------------------------------------------------------------------
// Main seed scanning logic (structures + biome checks)
bool scanSeed(uint64_t seed)
{
    bool hasAnyRequirements = false;
    bool allRequirementsMet = true;

    // Overworld generator
    Generator g;
    setupGenerator(&g, MC_1_21, 0);
    applySeed(&g, DIM_OVERWORLD, seed);

    // Determine bounding box around spawn or custom coords
    Pos spawn = {0,0};
    if (useSpawn) {
        spawn = getSpawn(&g);
    }
    int x0 = (useSpawn ? (spawn.x - searchRadius) : (customX - searchRadius));
    int z0 = (useSpawn ? (spawn.z - searchRadius) : (customZ - searchRadius));
    int x1 = x0 + (searchRadius * 2);
    int z1 = z0 + (searchRadius * 2);

    // Nether & End
    Generator ng, eg;
    setupGenerator(&ng, MC_1_21, 0);
    applySeed(&ng, DIM_NETHER, seed);
    setupGenerator(&eg, MC_1_21, 0);
    applySeed(&eg, DIM_END, seed);

    SurfaceNoise sn, esn;
    initSurfaceNoise(&sn, DIM_OVERWORLD, seed);
    initSurfaceNoise(&esn, DIM_END, seed);

    // 1) Structure cluster scanning
    if (clusterReq.enabled) {
        hasAnyRequirements = true;
        // Gather all structures for clusterReq.structureTypes within bounding box
        int capacity = 128;
        int clusterCount = 0;
        StructurePos *clusterPositions = malloc(capacity * sizeof(StructurePos));
        if (!clusterPositions) { perror("malloc"); exit(1); }

        for (int i = 0; i < clusterReq.count; i++) {
            int stype = clusterReq.structureTypes[i];
            StructureConfig sconf;
            if (!getStructureConfig(stype, MC_1_21, &sconf)) {
                // not valid in this version
                continue;
            }
            Generator *curr_gen = &g;
            SurfaceNoise *curr_sn = &sn;
            if (sconf.dim == DIM_NETHER) {
                curr_gen = &ng;
                curr_sn  = NULL; // if nether doesn't need surface noise
            }
            else if (sconf.dim == DIM_END) {
                curr_gen = &eg;
                curr_sn  = &esn;
            }

            double blocksPerRegion = sconf.regionSize * 16.0;
            int rx0 = (int)floor(x0 / blocksPerRegion);
            int rz0 = (int)floor(z0 / blocksPerRegion);
            int rx1 = (int)ceil(x1 / blocksPerRegion);
            int rz1 = (int)ceil(z1 / blocksPerRegion);

            // For each region, see if the structure can spawn
            for (int rz = rz0; rz <= rz1; rz++) {
                for (int rx = rx0; rx <= rx1; rx++) {
                    Pos pos;
                    if (!getStructurePos(stype, MC_1_21, seed, rx, rz, &pos))
                        continue;
                    if (pos.x < x0 || pos.x > x1 || pos.z < z0 || pos.z > z1)
                        continue;
                    if (!isViableStructurePos(stype, curr_gen, pos.x, pos.z, 0))
                        continue;

                    if (clusterCount == capacity) {
                        capacity *= 2;
                        clusterPositions = realloc(clusterPositions, capacity*sizeof(StructurePos));
                        if (!clusterPositions) { perror("realloc"); exit(1); }
                    }
                    clusterPositions[clusterCount].structureType = stype;
                    clusterPositions[clusterCount].x = pos.x;
                    clusterPositions[clusterCount].z = pos.z;
                    clusterCount++;
                }
            }
        }

        bool atLeastOneValidCluster = false;
        if (clusterCount >= 2) {
            // union-find by distance
            int *parent = malloc(clusterCount*sizeof(int));
            for (int i = 0; i < clusterCount; i++)
                parent[i] = i;

            for (int i = 0; i < clusterCount; i++) {
                for (int j = i+1; j < clusterCount; j++) {
                    int dx = clusterPositions[i].x - clusterPositions[j].x;
                    int dz = clusterPositions[i].z - clusterPositions[j].z;
                    // use squared distance compare
                    if (dx*dx + dz*dz <= clusterReq.clusterDistance * clusterReq.clusterDistance) {
                        unionSets(parent, i, j);
                    }
                }
            }

            bool *processed = calloc(clusterCount, sizeof(bool));
            for (int i = 0; i < clusterCount; i++) {
                int root = findSet(parent, i);
                if (processed[root]) continue;
                processed[root] = true;

                // gather all members in this cluster
                int *indices = malloc(16*sizeof(int));
                int indicesCap = 16;
                int groupSize = 0;

                for (int j = 0; j < clusterCount; j++) {
                    if (findSet(parent, j) == root) {
                        if (groupSize == indicesCap) {
                            indicesCap *= 2;
                            indices = realloc(indices, indicesCap*sizeof(int));
                        }
                        indices[groupSize++] = j;
                    }
                }

                // if cluster too small, skip
                if (groupSize < clusterReq.minClusterSize) {
                    free(indices);
                    continue;
                }

                // sort the structure types in ascending order
                int *groupTypes = malloc(groupSize*sizeof(int));
                for (int n = 0; n < groupSize; n++) {
                    groupTypes[n] = clusterPositions[indices[n]].structureType;
                }
                qsort(groupTypes, groupSize, sizeof(int), compareInts);

                // check for invalid combination
                bool invalid = isInvalidClusterDynamic(groupTypes, groupSize);
                if (invalid) {
                    printf("Skipping invalid cluster at seed %llu (contains an invalid combination)\n",
                           (unsigned long long) seed);
                    free(groupTypes);
                    free(indices);
                    continue;
                }

                // this cluster is valid
                atLeastOneValidCluster = true;

                // (Optional) print cluster info
                printf("== Seed %llu: Found cluster of size %d ==\n",
                       (unsigned long long) seed, groupSize);
                for (int n = 0; n < groupSize; n++) {
                    int idx = indices[n];
                    printf("   Type %d at (%d, %d)\n",
                           clusterPositions[idx].structureType,
                           clusterPositions[idx].x,
                           clusterPositions[idx].z);
                }
                printf("\n");
                //seedsFound++; // increment your global found count

                free(groupTypes);
                free(indices);
            }
            free(processed);
            free(parent);
        }
        free(clusterPositions);

        if (!atLeastOneValidCluster) {
            allRequirementsMet = false;
        }
    } // end if clusterReq.enabled

    // 2) Biome requirements
    if (biomeSearch.requiredCount > 0 || biomeSearch.clusterCount > 0) {
        hasAnyRequirements = true;
        if (!scanBiomes(&g, x0, z0, x1, z1, &biomeSearch)) {
            allRequirementsMet = false;
        }
    }

    // 3) Additional structure requirements array
    if (NUM_STRUCTURE_REQUIREMENTS > 0) {
        hasAnyRequirements = true;
        for (int rIndex = 0; rIndex < NUM_STRUCTURE_REQUIREMENTS; rIndex++) {
            StructureRequirement req = structureRequirements[rIndex];
            int foundCount = 0;

            StructureConfig sconf;
            if (!getStructureConfig(req.structureType, MC_1_21, &sconf)) {
                continue;
            }
            Generator *curr_gen = &g;
            SurfaceNoise *curr_sn = &sn;
            if (sconf.dim == DIM_NETHER) {
                curr_gen = &ng;
            }
            else if (sconf.dim == DIM_END) {
                curr_gen = &eg;
                curr_sn = &esn;
            }

            double blocksPerRegion = sconf.regionSize * 16.0;
            int rx0 = (int)floor(x0 / blocksPerRegion);
            int rz0 = (int)floor(z0 / blocksPerRegion);
            int rx1 = (int)ceil(x1 / blocksPerRegion);
            int rz1 = (int)ceil(z1 / blocksPerRegion);

            for (int rz = rz0; rz <= rz1; rz++) {
                for (int rx = rx0; rx <= rx1; rx++) {
                    Pos pos;
                    if (!getStructurePos(req.structureType, MC_1_21, seed, rx, rz, &pos))
                        continue;
                    if (pos.x < x0 || pos.x > x1 || pos.z < z0 || pos.z > z1)
                        continue;
                    if (!isViableStructurePos(req.structureType, curr_gen, pos.x, pos.z, 0))
                        continue;

                    // If we need to confirm biome is correct
                    int biome_id = getBiomeAt(curr_gen, 4, pos.x >> 2, pos.z >> 2, 320 >> 2);
                    if (biome_id == -1) {
                        // fallback
                        float heightArr[256];
                        int w = 16, h = 16;
                        Range r_range = {4, pos.x >> 2, pos.z >> 2, w, h, 320 >> 2, 1};
                        mapApproxHeight(heightArr, NULL, curr_gen, curr_sn, r_range.x, r_range.z, w, h);
                        int lx = pos.x & 15;
                        int lz = pos.z & 15;
                        int surface_y = (int)heightArr[lz*w + lx];
                        biome_id = getBiomeAt(curr_gen, 4, pos.x >> 2, surface_y >> 2, pos.z >> 2);
                    }

                    if (req.requiredBiome != -1 && biome_id != req.requiredBiome)
                        continue;

                    if (req.requiredBiome != -1 &&
                        (req.minBiomeSize != -1 || req.maxBiomeSize != -1))
                    {
                        int patchSize = getBiomePatchSize(curr_gen, pos.x, pos.z, biome_id);
                        if ((req.minBiomeSize != -1 && patchSize < req.minBiomeSize) ||
                            (req.maxBiomeSize != -1 && patchSize > req.maxBiomeSize))
                            continue;
                    }

                    foundCount++;
                    printf("Seed %llu: Found structure %d at (%d, %d) in biome %s\n",
                           (unsigned long long)seed, req.structureType, pos.x, pos.z,
                           getBiomeName(biome_id));
                    //seedsFound++;
                }
            }

            // If no structures found that match this requirement => fail
            if (foundCount < req.minCount) {
                allRequirementsMet = false;
            }
        }
    }

    if (!hasAnyRequirements) {
        printf("Warning: No requirements set, skipping validation for seed %llu.\n",
               (unsigned long long) seed);
        return false;
    }

    if (allRequirementsMet) {
        printf("Valid seed found: %llu\n", (unsigned long long) seed);
        return true;
    }
    return false;
}

// -----------------------------------------------------------------------------
// Thread function: scans seeds in [currentSeed..end_seed]
void *scanTask(void *arg)
{
    uint64_t *endSeedPtr = (uint64_t*) arg;

    while (true) {
        pthread_mutex_lock(&seedMutex);

        // Remove or comment out these lines so the loop never breaks prematurely:
        // if (foundValidSeed || currentSeed > *endSeedPtr) {
        //     pthread_mutex_unlock(&seedMutex);
        //     break;
        // }

        uint64_t seed = currentSeed;
        currentSeed++;
        pthread_mutex_unlock(&seedMutex);

        // Now scan that seed
        if (scanSeed(seed)) {
            pthread_mutex_lock(&seedMutex);
            // Optionally re-add increment if you want to track how many valid seeds were found:
            // seedsFound++;

            validSeed = seed;
            // Remove or comment out the forced exit so it won't terminate all threads:
            // if (seedsFound >= MAX_SEEDS_TO_FIND) {
            //     foundValidSeed = true;
            //     pthread_mutex_unlock(&seedMutex);
            //     exit(0);
            // }

            pthread_mutex_unlock(&seedMutex);

            // Likewise, do not break here unless you want each thread to stop after first success:
            // break;
        }
    }
    return NULL;
}

// -----------------------------------------------------------------------------
// Main function: sets up multi-thread scanning
int main()
{
    printf("Checking cubiomes library...\n");

    // For demonstration, let's scan 1000 seeds starting at starting_seed
    uint64_t scanCount = 10000;
    end_seed = starting_seed + scanCount - 1ULL;

    // Set number of threads:
    tasksCount = 100; // e.g. 2 threads. You can set this to whatever you like.

    // Initialize the Biome Requirement structures
    reqGroup.biomeIds      = reqGroup0;
    reqGroup.biomeCount    = sizeof(reqGroup0)/sizeof(reqGroup0[0]);
    reqGroup.sizeConfigs   = reqSizeConfigs;
    reqGroup.configCount   = sizeof(reqSizeConfigs)/sizeof(reqSizeConfigs[0]);
    reqGroup.logCenters    = 1;
    requiredBiomes[0]      = reqGroup;

    // Initialize global scanning state
    currentSeed     = starting_seed;
    foundValidSeed  = false;
    seedsFound      = 0;

    // Create threads
    pthread_t *threads = malloc(tasksCount * sizeof(pthread_t));
    for (int i = 0; i < tasksCount; i++) {
        // pass address of end_seed so each thread knows the limit
        pthread_create(&threads[i], NULL, scanTask, &end_seed);
    }

    // Wait for all threads to finish
    for (int i = 0; i < tasksCount; i++) {
        pthread_join(threads[i], NULL);
    }
    free(threads);

    // Final result
    if (foundValidSeed) {
        printf("== Found at least one valid seed (e.g., %llu) ==\n",
               (unsigned long long) validSeed);
    }
    else {
        printf("Finished searching, no valid seeds found in [%llu..%llu].\n",
               (unsigned long long)starting_seed, (unsigned long long)end_seed);
    }

    return 0;
}
