# **Advanced Minecraft Seed Finder**  

## **Website:** [mc-seed-finder.replit.app](https://mc-seed-finder.replit.app/)  

### *(If website is down, you can use the seed finder by following the Installation)*

## **Overview**  
This tool automates finding Minecraft seeds based on specific biomes, structures, and terrain features. Supports Java and Bedrock editions (1.18 to 1.21.6)

## **Features**  
- Find structures with biome constraints, height ranges, and minimum amount  
- Find any biomes and any biome combinations
- Locate clustered structures (structures next to each other)
- Set biome size requirements

## **Usage**  
After running, set and paste the parameters for seed scanning, and then press Control+D enter the parameters to start the scan.

### **Working Example:**
```
===== Scanning options =====
Starting seed: 0
Search radius: 1000
Use spawn: true
Custom z: 0
Custom x: 0

===== Required structures =====

===== Structure Clusters =====
Enabled: false
Valid biomes: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 16, 18
Min cluster distance = 32,
Min cluster size  = 2

===== Invalid Structure Clusters =====
1. 16, 11
2. 16, 6

===== Required Biomes =====
1. 185 (min size: -1, max size: -1)

===== Clustered Biomes =====
```
### **Another Example:**
```
===== Scanning options =====
Starting seed: 0
Search radius: 1000
Use spawn: false
Custom z: 0
Custom x: 0

===== Required structures =====
1. 9 (min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1)
2. 5 (min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1)

===== Structure Clusters =====
Enabled: true
Valid biomes: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 16, 18
Min cluster distance = 32,
Min cluster size  = 2

===== Invalid Structure Clusters =====
1. 16, 11
2. 16, 6

===== Required Biomes =====
1. 1 (min size: 50, max size: -1)
2. 185 (min size: -1, max size: -1)
3. 4 (min size: -1, max size: -1)

===== Clustered Biomes =====
1. 185, 1 (min size: -1, max size: -1)
2. 4, 1 (min size: 100, max size: -1)
2. 185, 4 (min size: -1, max size: -1)
```

## **Installation**  

### **Clone the Repository**  
```bash
git clone https://github.com/Dalle2fan/cubiomes.git
```
```bash
cd cubiomes
```
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
```bash
brew install ruby
brew install cmake
```
### **Run Application**  
```bash
chmod +x build.sh
sudo ./build.sh
```

## **Uninstallation**
```bash
sudo chmod -R u+rwx cubiomes
sudo rm -rf cubiomes
```
### **Biome and Structures ID references:**
```
// Structure ID Reference
//  0  - Feature
//  1  - Desert_Pyramid
//  2  - Jungle_Temple
//  3  - Swamp_Hut
//  4  - Igloo
//  5  - Village
//  6  - Ocean_Ruin
//  7  - Shipwreck
//  8  - Monument
//  9  - Mansion
// 10  - Outpost
// 11  - Ruined_Portal
// 12  - Ruined_Portal_N
// 13  - Ancient_City
// 14  - Treasure
// 15  - Mineshaft
// 16  - Desert_Well
// 17  - Geode
// 18  - Trail_Ruins
// 19  - Trial_Chambers

// Biome ID Reference
//  0   - Ocean
//  1   - Plains
//  2   - Desert
//  3   - Windswept Hills
//  4   - Forest
//  5   - Taiga
//  6   - Swamp
//  7   - River
// 10   - Frozen Ocean
// 11   - Frozen River
// 12   - Snowy Plains
// 13   - Snowy Mountains
// 14   - Mushroom Fields
// 15   - Mushroom Fields Shore
// 16   - Beach
// 17   - Desert Hills
// 18   - Windswept Forest
// 19   - Taiga Hills
// 20   - Mountain Edge
// 21   - Jungle
// 22   - Jungle Hills
// 23   - Sparse Jungle
// 24   - Deep Ocean
// 25   - Stony Shore
// 26   - Snowy Beach
// 27   - Birch Forest
// 28   - Birch Forest Hills
// 29   - Dark Forest
// 30   - Snowy Taiga
// 31   - Snowy Taiga Hills
// 32   - Old Growth Pine Taiga
// 33   - Giant Tree Taiga Hills
// 34   - Wooded Mountains
// 35   - Savanna
// 36   - Savanna Plateau
// 37   - Badlands
// 38   - Wooded Badlands
// 39   - Badlands Plateau
// 44   - Warm Ocean
// 45   - Lukewarm Ocean
// 46   - Cold Ocean
// 47   - Deep Warm Ocean
// 48   - Deep Lukewarm Ocean
// 49   - Deep Cold Ocean
// 50   - Deep Frozen Ocean
// 129  - Sunflower Plains
// 130  - Desert Lakes
// 131  - Windswept Gravelly Hills
// 132  - Flower Forest
// 133  - Taiga Mountains
// 134  - Swamp Hills
// 140  - Ice Spikes
// 149  - Modified Jungle
// 151  - Modified Jungle Edge
// 155  - Old Growth Birch Forest
// 156  - Tall Birch Hills
// 157  - Dark Forest Hills
// 158  - Snowy Taiga Mountains
// 160  - Old Growth Spruce Taiga
// 161  - Giant Spruce Taiga Hills
// 162  - Gravelly Mountains+
// 163  - Windswept Savanna
// 164  - Shattered Savanna Plateau
// 165  - Eroded Badlands
// 166  - Modified Wooded Badlands Plateau
// 167  - Modified Badlands Plateau
// 168  - Bamboo Jungle
// 169  - Bamboo Jungle Hills
// 170  - Underground
// 171  - Underground
// 172  - Underground
// 174  - Dripstone Caves
// 175  - Lush Caves
// 177  - Meadow
// 178  - Grove
// 179  - Snowy Slopes
// 180  - Frozen Peaks
// 181  - Jagged Peaks
// 182  - Stony Peaks
// 183  - Deep Dark
// 184  - Mangrove Swamp
// 185  - Cherry Grove
// 186  - Pale Garden
```
