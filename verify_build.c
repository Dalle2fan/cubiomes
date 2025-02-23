// verify_build.c
#include "cubiomes/finders.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>
#include <string.h> // For memcpy, strcmp, etc.

// -----------------------------------------------------------------------------
// Global configuration / defaults
#define MAX_SEEDS_TO_FIND 1 
int seedsFound = 0;          // Tracks how many seeds have been found so far

// Range for seed scanning:
uint64_t starting_seed = 1;
uint64_t end_seed      = 0ULL;  // Will be set in main (starting_seed + some offset)

/* Uses this:
===== Scanning options =====
Starting seed: 0
Search radius: 1000
Use spawn: false
Custom z: 0
Custom x: 0
    */
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
// Biome ID -> Name helper
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
                // slight optimization: break if we’re well past the region
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
    int structureType;   // e.g. 5 -> Village, etc. (Cubiomes structure ID)
    int minCount;
    int minHeight;
    int maxHeight;
    int requiredBiome;   // -1 => skip biome check
    int minBiomeSize;    // -1 => no minimum
    int maxBiomeSize;    // -1 => no maximum
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
    int logCenters; // whether to print out patch centers
} BiomeRequirement;

// For "clustered" biomes
typedef struct {
    int *biomeIds;
    int biomeCount;
    int minSize;
    int maxSize;
    int logCenters; // whether to print out cluster centers
} BiomeCluster;

// Top-level container for biome searches
typedef struct {
    BiomeRequirement *required;
    int requiredCount;
    BiomeCluster *clusters;
    int clusterCount;
} BiomeSearch;

// We will dynamically build these from the parameter file
BiomeRequirement *g_requiredBiomes   = NULL; 
int g_requiredBiomesCount            = 0;

BiomeCluster    *g_biomeClusters     = NULL;
int g_biomeClustersCount             = 0;

BiomeSearch biomeSearch = {
    .required      = NULL,
    .requiredCount = 0,
    .clusters      = NULL,
    .clusterCount  = 0
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

// We will fill this from the parameter file
ClusterRequirement clusterReq = {
    .enabled         = false,
    .clusterDistance = 32,
    .structureTypes  = NULL,
    .count           = 0,
    .minClusterSize  = 2
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

// We’ll store them dynamically as well
InvalidCombination *invalidCombinations = NULL;
int numInvalidCombinations = 0;

// -----------------------------------------------------------------------------
// Our array of required structures (dynamically built)
StructureRequirement *structureRequirements = NULL;
int NUM_STRUCTURE_REQUIREMENTS = 0;

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
    for (int ic = 0; ic < numInvalidCombinations; ic++) {
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
    if (bs->requiredCount == 0 && bs->clusterCount == 0)
    {
        // If no biome requirements are present, we consider it "passing" or you can handle otherwise
        return true;
    }
    bool success = true; 
    int step = 4;

    // 1) Required biome patches
    if (bs->requiredCount > 0) 
    {
        bool allRequirementsFound = true;

        // We'll require that at least one patch for *each* BiomeRequirement is found
        // If your logic is different (like you need them all simultaneously?), you can tweak.
        for (int i = 0; i < bs->requiredCount; i++) {
            BiomeRequirement *req = &bs->required[i];
            bool foundPatchForThisReq = false;

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
                allRequirementsFound = false;
                continue;
            }

            // union-find
            int *parent = malloc(count * sizeof(int));
            if (!parent) { perror("malloc"); exit(1); }
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
            if (!processed) { perror("calloc"); exit(1); }

            // For each root in union-find, check if it meets size constraints
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
                    foundPatchForThisReq = true;
                    if (req->logCenters) {
                        printf("Required biome patch (reqIndex=%d, biome=%s): center (%.1f,%.1f), cells=%d\n",
                               i, getBiomeName(theBiome), centerX, centerZ, compCount);
                    }
                }
            }
            free(processed);
            free(parent);
            free(positions);

            if (!foundPatchForThisReq) {
                allRequirementsFound = false;
            }
        } // end for each required group

        if (!allRequirementsFound) {
            success = false;
        }
    }

    // 2) Clustered biomes
    if (bs->clusterCount > 0) 
    {
        // For each cluster definition, check if we can find it
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
                // You might decide that not finding any means failure or not,
                // depending on your needs. Here, we just continue.
                continue;
            }

            // union-find them by adjacency
            int *parent = malloc(count * sizeof(int));
            if (!parent) { perror("malloc"); exit(1); }
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
            if (!processed) { perror("calloc"); exit(1); }

            for (int c = 0; c < count; c++) {
                int root = findSet(parent, c);
                if (processed[root]) continue;
                processed[root] = true;

                double sumX = 0, sumZ = 0;
                int compCount = 0;
                // Track which biome IDs appear
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

                // If you specifically want pairs or a certain pattern, do that check here.
                // For demonstration, we just check min/max size
                bool sizeOk = true;
                if (cl->minSize > -1 && compCount < cl->minSize) sizeOk = false;
                if (cl->maxSize > -1 && compCount > cl->maxSize) sizeOk = false;

                if (sizeOk && cl->logCenters) {
                    double cx = sumX / compCount;
                    double cz = sumZ / compCount;
                    printf("Clustered biome group %d: (distinct biome count=%d), center(%.1f,%.1f), count=%d\n",
                           i+1, dIdx, cx, cz, compCount);
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
        // Gather all structures (of the specified cluster types) within bounding box
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
                curr_sn  = NULL; // Nether typically doesn’t need it here
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

                    // Save it
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
            if (!parent) { perror("malloc"); exit(1); }
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
            if (!processed) { perror("calloc"); exit(1); }

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

                    // Check biome only if requiredBiome != -1
                    int biome_id = getBiomeAt(curr_gen, 4, pos.x >> 2, pos.z >> 2, 320 >> 2);
                    if (biome_id == -1) {
                        // fallback with approximate height
                        float heightArr[256];
                        int w = 16, h = 16;
                        Range r_range = {4, pos.x >> 2, pos.z >> 2, w, h, 320 >> 2, 1};
                        mapApproxHeight(heightArr, NULL, curr_gen, curr_sn, r_range.x, r_range.z, w, h);
                        int lx = pos.x & 15;
                        int lz = pos.z & 15;
                        int surface_y = (int)heightArr[lz*w + lx];
                        biome_id = getBiomeAt(curr_gen, 4, surface_y >> 2, surface_y >> 2, pos.z >> 2);
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

                    // Height check
                    // (If your code uses terrain-based Y, you'll need to sample it. 
                    //  Here we only do a simple pass, so let's skip. 
                    //  But you can do something like:
                    // int terrainY = estimateSurfaceHeight(...);
                    // if (terrainY < req.minHeight || terrainY > req.maxHeight) continue;
                    //)

                    foundCount++;
                    printf("Seed %llu: Found structure %d at (%d, %d) in biome %s\n",
                           (unsigned long long)seed, req.structureType, pos.x, pos.z,
                           getBiomeName(biome_id));
                }
            }

            // If no structures found that match this requirement => fail
            if (foundCount < req.minCount) {
                allRequirementsMet = false;
            }
        }
    }

    if (!hasAnyRequirements) {
        // If absolutely no requirements are set, you can decide to skip or treat it as valid.
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
        // If you want to stop scanning after reaching end_seed, do:
        if (currentSeed > *endSeedPtr)
        {
            pthread_mutex_unlock(&seedMutex);
            break;
        }
        uint64_t seed = currentSeed;
        currentSeed++;
        pthread_mutex_unlock(&seedMutex);

        // Now scan that seed
        if (scanSeed(seed)) {
            pthread_mutex_lock(&seedMutex);
            validSeed = seed;
            seedsFound++;
            // If you want to stop once we find enough seeds, do something like:
            if (seedsFound >= MAX_SEEDS_TO_FIND) {
                foundValidSeed = true;
                pthread_mutex_unlock(&seedMutex);
                // optionally exit the entire process:
                // exit(0);
                // or break to let just this thread stop scanning:
                break;
            }
            pthread_mutex_unlock(&seedMutex);
        }
    }
    return NULL;
}

// -----------------------------------------------------------------------------
// Simple helper to skip leading/trailing spaces
static void trim(char *str)
{
    // left trim
    char *p = str;
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    if (p != str) memmove(str, p, strlen(p)+1);

    // right trim
    int len = (int)strlen(str);
    while (len > 0 && (str[len-1] == ' ' || str[len-1] == '\t' 
                       || str[len-1] == '\r' || str[len-1] == '\n')) {
        str[len-1] = '\0';
        len--;
    }
}

// -----------------------------------------------------------------------------
// Parsing function to read the config from lines (file or otherwise).
// This is a simplistic parser that looks for known headings and then key lines.
void parseConfigLine(const char *section, char *line)
{
    // Example lines:
    // "Starting seed: 0"
    // "Search radius: 1000"
    // "Use spawn: false"
    // ...
    if (strcmp(section, "===== Scanning options =====") == 0) 
    {
        // We'll parse scanning options
        if (strstr(line, "Starting seed:") == line) {
            sscanf(line, "Starting seed: %" SCNu64, &starting_seed);
        }
        else if (strstr(line, "Search radius:") == line) {
            sscanf(line, "Search radius: %d", &searchRadius);
        }
        else if (strstr(line, "Use spawn:") == line) {
            char val[16];
            if (sscanf(line, "Use spawn: %15s", val) == 1) {
                // interpret "true"/"1" as 1, else 0
                if (strcmp(val, "true") == 0 || strcmp(val, "1") == 0) 
                    useSpawn = 1;
                else
                    useSpawn = 0;
            }
        }
        else if (strstr(line, "Custom x:") == line) {
            sscanf(line, "Custom x: %d", &customX);
        }
        else if (strstr(line, "Custom z:") == line) {
            sscanf(line, "Custom z: %d", &customZ);
        }
    }
    // etc. For other sections, do similarly.
}

// A few global variables to track which section we're parsing
static char currentSection[128] = {0};

void parseParameterLine(char *line)
{
    trim(line);
    if (line[0] == '\0') {
        return; // skip empty
    }

    // Check if it's a section heading line
    if (strstr(line, "=====") == line) {
        // This means we've reached a new section
        strcpy(currentSection, line);
        return;
    }

    // Otherwise, parse based on the currentSection
    if (strcmp(currentSection, "===== Scanning options =====") == 0) 
    {
        parseConfigLine(currentSection, line);
    }
    else if (strcmp(currentSection, "===== Required structures =====") == 0) 
    {
        // Lines look like:
        // 1. Village (min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1)
        // We just parse them
        // We might do something like:
        int idx, minCount, minH, maxH, biome, minSz, maxSz;
        char nameBuf[64];
        // Because the line might have the structure name in parentheses, 
        // we can do a rough parse. Adjust as needed for your actual structure IDs.
        // Example line:
        // "2. Mansion (min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1)"
        // We'll pretend we know how to map "Mansion" -> Cubiomes ID 6 or something.
        // For demonstration, let's do a trivial approach:

        // Find structure name: everything up to '(' minus the "x. " prefix
        // This requires a robust approach or a known set of structure names -> IDs.

        // For brevity, assume we can parse:
        //   <idx>. <Name> (min amount: ...
        // Then we look inside parentheses.

        // This is an example. Real code might have a lookup table from "Village" to 5, "Mansion" to 6, etc.
        // Let's do a small switch as a demo:
        // We'll define a helper function or inline code.

        // Quick parse (this is *very* naive):
        // Step 1: parse out index and name:

        char *openParen = strchr(line, '(');
        if (!openParen) return; 
        // parse "x. <Name>" part
        // e.g. "1. Village "
        char prefix[128];
        strncpy(prefix, line, openParen - line);
        prefix[openParen-line] = '\0';
        trim(prefix);

        // prefix might be "1. Village"
        sscanf(prefix, "%d. %63[^\n]", &idx, nameBuf);
        trim(nameBuf);

        // Step 2: parse parentheses
        // e.g. "min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1)"
        char parenPart[256];
        strcpy(parenPart, openParen+1);
        // remove trailing ')'
        char *endParen = strrchr(parenPart, ')');
        if (endParen) *endParen = '\0';
        // Now parse the known fields:
        // "min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1"
        sscanf(parenPart, "min amount: %d, min height: %d, max height: %d, biome: %d, min size: %d, max size: %d",
                          &minCount, &minH, &maxH, &biome, &minSz, &maxSz);

        // Map the nameBuf to actual structure ID
        // Just as an example:
        int structureType = -1;
        if (strstr(nameBuf, "Village")) {
            // For cubiomes: 5 is Village (in the structure enum)
            structureType = 5; 
        }
        else if (strstr(nameBuf, "Mansion")) {
            // 6 is Mansion in the structure enum? (Double-check your actual ID.)
            structureType = 6;
        }
        else {
            // add more mappings as needed
        }

        // Append to structureRequirements array
        structureRequirements = realloc(structureRequirements, (NUM_STRUCTURE_REQUIREMENTS+1)*sizeof(StructureRequirement));
        structureRequirements[NUM_STRUCTURE_REQUIREMENTS].structureType = structureType;
        structureRequirements[NUM_STRUCTURE_REQUIREMENTS].minCount       = minCount;
        structureRequirements[NUM_STRUCTURE_REQUIREMENTS].minHeight      = minH;
        structureRequirements[NUM_STRUCTURE_REQUIREMENTS].maxHeight      = maxH;
        structureRequirements[NUM_STRUCTURE_REQUIREMENTS].requiredBiome  = biome;
        structureRequirements[NUM_STRUCTURE_REQUIREMENTS].minBiomeSize   = minSz;
        structureRequirements[NUM_STRUCTURE_REQUIREMENTS].maxBiomeSize   = maxSz;
        NUM_STRUCTURE_REQUIREMENTS++;
    }
    else if (strcmp(currentSection, "===== Structure Clusters =====") == 0)
    {
        // Example lines:
        //   Enabled: true
        //   Valid biomes: 1, 2, 3, 4, ...
        //   Min cluster distance = 32,
        //   Min cluster size  = 2
        // You can parse them similarly
        if (strstr(line, "Enabled:") == line) {
            char val[16];
            if (sscanf(line, "Enabled: %15s", val) == 1) {
                if (strcmp(val, "true") == 0 || strcmp(val, "1") == 0) 
                    clusterReq.enabled = true;
                else 
                    clusterReq.enabled = false;
            }
        }
        else if (strstr(line, "Valid biomes:") == line) {
            // Actually, the example said "Valid biomes: 1,2,3..." 
            // but the code snippet you gave uses structure IDs, not biome IDs, 
            // so let's assume these are the structure *types* we want to cluster:
            // We'll parse them into clusterReq.structureTypes
            const char *p = strchr(line, ':');
            if (p) {
                p++; // skip colon
                // parse comma-separated
                char list[256];
                strcpy(list, p);
                trim(list);
                // tokenize
                char *tok = strtok(list, ",");
                while (tok) {
                    int stype = atoi(tok);
                    clusterReq.structureTypes = realloc(clusterReq.structureTypes, (clusterReq.count+1)*sizeof(int));
                    clusterReq.structureTypes[clusterReq.count] = stype;
                    clusterReq.count++;
                    tok = strtok(NULL, ",");
                }
            }
        }
        else if (strstr(line, "Min cluster distance") == line) {
            int dist;
            if (sscanf(line, "Min cluster distance = %d", &dist) == 1) {
                clusterReq.clusterDistance = dist;
            }
        }
        else if (strstr(line, "Min cluster size") == line) {
            int sz;
            if (sscanf(line, "Min cluster size = %d", &sz) == 1) {
                clusterReq.minClusterSize = sz;
            }
        }
    }
    else if (strcmp(currentSection, "===== Invalid Structure Clusters =====") == 0)
    {
        // Example lines:
        //   1. 16, 11
        //   2. 16, 6
        // parse each line as a list of structure IDs
        // We'll store them in invalidCombinations
        char *dot = strchr(line, '.');
        if (!dot) return;
        dot++; // skip the dot
        trim(dot);
        // dot might be "16, 11"
        // parse them
        int tmpArr[64];
        int tmpCount = 0;
        char *tok = strtok(dot, ",");
        while (tok) {
            tmpArr[tmpCount++] = atoi(tok);
            tok = strtok(NULL, ",");
        }
        // store in invalidCombinations
        invalidCombinations = realloc(invalidCombinations, (numInvalidCombinations+1)*sizeof(InvalidCombination));
        invalidCombinations[numInvalidCombinations].count = tmpCount;
        invalidCombinations[numInvalidCombinations].types = malloc(tmpCount*sizeof(int));
        for (int i = 0; i < tmpCount; i++) {
            invalidCombinations[numInvalidCombinations].types[i] = tmpArr[i];
        }
        numInvalidCombinations++;
    }
    else if (strcmp(currentSection, "===== Required Biomes =====") == 0)
    {
        // Example lines:
        // 1. 1 (min size: 50, max size: -1)
        // 2. 185 (min size: -1, max size: -1)
        // ...
        // We’ll dynamically build g_requiredBiomes as an array of BiomeRequirement,
        // each with a single biome or potentially more, depending on your structure.

        // But from your example, each line is 1 biome + min/max size.
        // We might store them separately and then finalize them after parsing
        static BiomeSizeConfig *tmpConfigs = NULL;
        static int tmpConfigsCount = 0;

        int idx, biomeId, minSz, maxSz;
        // parse line e.g.: "1. 1 (min size: 50, max size: -1)"
        char *paren = strchr(line, '(');
        if (!paren) return;
        // parse "x. y"
        // e.g. "1. 1"
        int consumed = 0;
        sscanf(line, "%d. %d%n", &idx, &biomeId, &consumed);

        // parse parentheses
        char parenPart[256];
        strcpy(parenPart, paren+1);
        char *endParen = strrchr(parenPart, ')');
        if (endParen) *endParen = '\0';
        // "min size: 50, max size: -1"
        sscanf(parenPart, "min size: %d, max size: %d", &minSz, &maxSz);

        // We'll store them in a temporary array. 
        tmpConfigs = realloc(tmpConfigs, (tmpConfigsCount+1)*sizeof(BiomeSizeConfig));
        tmpConfigs[tmpConfigsCount].biomeId = biomeId;
        tmpConfigs[tmpConfigsCount].minSize = minSz;
        tmpConfigs[tmpConfigsCount].maxSize = maxSz;
        tmpConfigsCount++;

        // For demonstration, let's keep them all in a single "required group".
        // Or you might want each line to be a separate requirement. 
        // For simplicity, let's do *one* BiomeRequirement that includes multiple biome IDs.

        // If your logic requires each line to be its own BiomeRequirement, you'd do so here.

        // We'll create or expand a single BiomeRequirement with all these. 
        if (g_requiredBiomesCount == 0) {
            g_requiredBiomes = malloc(sizeof(BiomeRequirement));
            g_requiredBiomesCount = 1;
            g_requiredBiomes[0].biomeIds     = NULL;
            g_requiredBiomes[0].biomeCount   = 0;
            g_requiredBiomes[0].sizeConfigs  = NULL;
            g_requiredBiomes[0].configCount  = 0;
            g_requiredBiomes[0].logCenters   = 1; // default
        }
        BiomeRequirement *br = &g_requiredBiomes[0];
        // add this biomeId to br->biomeIds
        br->biomeIds = realloc(br->biomeIds, (br->biomeCount+1)*sizeof(int));
        br->biomeIds[br->biomeCount] = biomeId;
        br->biomeCount++;

        // Also store the size constraints. 
        // In a more robust design, you'd match them by index or by biome ID. 
        // For now, let's just store them all in a single array:
        br->sizeConfigs = realloc(br->sizeConfigs, tmpConfigsCount*sizeof(BiomeSizeConfig));
        memcpy(br->sizeConfigs, tmpConfigs, tmpConfigsCount*sizeof(BiomeSizeConfig));
        br->configCount = tmpConfigsCount;
    }
    else if (strcmp(currentSection, "===== Clustered Biomes =====") == 0)
    {
        // Example lines:
        // 1. 185, 1 (min size: -1, max size: -1)
        // 2. 4, 1 (min size: 100, max size: -1)
        // ...
        // We'll parse similarly. Suppose each line is one cluster definition.
        // Then we store it into g_biomeClusters.

        int idx, minSz, maxSz;
        int b1, b2;
        // naive parse: "1. 185, 1 (min size: -1, max size: -1)"
        // or "2. 4, 1 (min size: 100, max size: -1)"
        char *paren = strchr(line, '(');
        if (!paren) return;

        char prefix[128];
        strncpy(prefix, line, paren - line);
        prefix[paren - line] = '\0';
        trim(prefix);
        // e.g. "1. 185, 1"
        sscanf(prefix, "%d. %d, %d", &idx, &b1, &b2);

        // parse parentheses
        char parenPart[128];
        strcpy(parenPart, paren+1);
        char *endParen = strrchr(parenPart, ')');
        if (endParen) *endParen = '\0';
        trim(parenPart);
        // "min size: -1, max size: -1"
        sscanf(parenPart, "min size: %d, max size: %d", &minSz, &maxSz);

        // Now store in g_biomeClusters
        g_biomeClusters = realloc(g_biomeClusters, (g_biomeClustersCount+1)*sizeof(BiomeCluster));
        BiomeCluster *bc = &g_biomeClusters[g_biomeClustersCount];
        g_biomeClustersCount++;

        bc->biomeIds = malloc(2*sizeof(int)); // we found 2 from the line
        bc->biomeIds[0] = b1;
        bc->biomeIds[1] = b2;
        bc->biomeCount  = 2;
        bc->minSize     = minSz;
        bc->maxSize     = maxSz;
        bc->logCenters  = 1; // or false, if you prefer
    }
}

// -----------------------------------------------------------------------------
// A function to read from a file or from stdin, building our config
void parseParameterStream(FILE *fp)
{
    char buf[512];
    while (fgets(buf, sizeof(buf), fp)) {
        parseParameterLine(buf);
    }
}

// -----------------------------------------------------------------------------
// MAIN
int main(int argc, char *argv[])
{
    printf("=== Parameter-Based Scanning ===\n");
    printf("Attach a .txt file or specify a path to the config, or paste config lines.\n");

    // 1) Try to read a file if given:
    FILE *fp = NULL;
    if (argc > 1) {
        fp = fopen(argv[1], "r");
        if (!fp) {
            fprintf(stderr, "Could not open file '%s'. Using stdin instead.\n", argv[1]);
        }
    }
    if (!fp) {
        printf("Please paste your config (end with Ctrl+D or Ctrl+Z on Windows):\n");
        fp = stdin;
    }

    // 2) Parse the config
    parseParameterStream(fp);
    if (fp != stdin) fclose(fp);

    // 3) Now we have our global variables from the config, set up the BiomeSearch, etc.
    biomeSearch.required      = g_requiredBiomes;
    biomeSearch.requiredCount = g_requiredBiomesCount;
    biomeSearch.clusters      = g_biomeClusters;
    biomeSearch.clusterCount  = g_biomeClustersCount;

    // 4) Decide how many seeds to scan
    // For demonstration, let's just scan 10k seeds from starting_seed
    uint64_t scanCount = 10000;
    end_seed = starting_seed + scanCount - 1ULL;

    // 5) Set number of threads (you might parse this from the config too)
    // For example, let's keep it at 2
    tasksCount = 1;

    // 6) Initialize global scanning state
    currentSeed     = starting_seed;
    foundValidSeed  = false;
    seedsFound      = 0;

    // 7) Spawn threads
    pthread_t *threads = malloc(tasksCount * sizeof(pthread_t));
    for (int i = 0; i < tasksCount; i++) {
        pthread_create(&threads[i], NULL, scanTask, &end_seed);
    }

    // 8) Wait for all threads to finish
    for (int i = 0; i < tasksCount; i++) {
        pthread_join(threads[i], NULL);
    }
    free(threads);

    // 9) Final result
    if (foundValidSeed) {
        printf("== Found at least one valid seed (e.g., %llu). Seeds found: %d ==\n",
               (unsigned long long) validSeed, seedsFound);
    }
    else {
        printf("Finished searching, no valid seeds found in [%llu..%llu].\n",
               (unsigned long long)starting_seed, (unsigned long long)end_seed);
    }

    // Clean up your dynamic allocations if you like
    // e.g. free(g_requiredBiomes), free(invalidCombinations[x].types), etc.

    return 0;
}
