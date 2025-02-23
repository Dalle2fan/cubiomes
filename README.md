# **Advanced Minecraft Seed Finder**  

## **Website:** [mc-seed-finder.replit.app](https://mc-seed-finder.replit.app/)  

### *(If website is down, you can use the seed finder by following the Installation)*

## **Overview**  
This tool automates finding Minecraft seeds based on specific biomes, structures, and terrain features. Supports Java and Bedrock editions (1.18 to 1.21.6)

## **Features**  
- Find any biomes including custom biomes (like islands, valleys, and elevated encircling terrain)
- Find structures with biome constraints, height ranges, and minimum counts  
- Locate clustered structures (structures next to each other) and biome combinations  
- Set biome size requirements for more precise world generation  

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
1. Village (min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1)
2. Mansion (min amount: 1, min height: -9999, max height: 9999, biome: -1, min size: -1, max size: -1)

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
