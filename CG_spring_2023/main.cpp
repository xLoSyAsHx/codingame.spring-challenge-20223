#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <queue>
#include <map>
#include <stack>
#include <type_traits>

#include <iomanip>
#include <memory>
#include <random>
#include <numeric>
#include <algorithm>

#define DBG_MSG_S(s) { std::cerr << s << std::endl; }

#define DBG_MSG_V(msg, v) { std::cerr << msg << " " << v << std::endl; }
#define DBG_MSG_V2(msg, v1, v2) { std::cerr << msg << " "#v1 << ' ' << v1 << " " << #v2 << ' ' << v2 << std::endl; }
#define DBG_MSG_V3(msg, v1, v2, v3) { std::cerr << msg << " " << #v1 << ' ' << v1 << " " << #v2 << ' ' << v2 << " " << #v3 << ' ' << v3 << std::endl; }

#define DBG_V(v) { std::cerr << ""#v" = " << v << std::endl; }
#define DBG_V2(v1, v2) { std::cerr << " " << #v1 << ' ' << v1 << "; " << #v2 << ' ' << v2 << std::endl; }
#define DBG_V3(v1, v2, v3) { std::cerr << " " << #v1 << ' ' << v1 << ";  " << #v2 << ' ' <<  v2 << "; " << #v3 << ' ' << v3 <<  std::endl; }
#define DBG_V4(v1, v2, v3, v4) { std::cerr << " " << #v1 << ' ' << v1 << ";  " << #v2 << ' ' <<  v2 << "; " << #v3 << ' ' << v3 << "; " << #v4 << ' ' << v4 <<  std::endl; }

#define DBG_MSG_ARR_V(v) { std::cerr << ""#v" = { "; for (auto& el : v) std::cerr << el << "; "; std::cerr << std::endl; }

enum MC_Type {
    MC_EMPTY = 0,
    MC_EGGS = 1,
    MC_CRYSTAL = 2,
};

struct MapCell {
    std::array<MapCell*, 6> neighbours;

    bool isValid = false;

    MC_Type type;
    int idx;
    int resources;

    int my_ants;
    int opp_ants;

    // debug

    // Global ptr to map
    static std::vector<MapCell*>& map;
};

struct GameInfo
{
    std::mt19937 random_device{ std::random_device{}() };

    const int MAX_CRYSTALS_SPOT_DISTANCE = 3;
    const int MAX_EGGS_SPOT_DISTANCE = 2;

    int currentTurn = -1;

    int init_crystals = 0;
    int cur_crystals = 0;
    int init_eggs = 0;
    int cur_eggs = 0;

    int my_ants = 0;
    int opp_ants = 0;

    int num_of_bases = 0;
    int number_of_cells = 0;
    std::vector<MapCell*> map = {};
    std::vector<MapCell*> crystals = {};
    std::vector<MapCell*> eggs = {};
    std::vector<MapCell*> my_bases = {};
    std::vector<MapCell*> opp_bases = {};

    void Initialize()
    {

        std::cin >> number_of_cells; std::cin.ignore();

        map.resize(number_of_cells);
        for (int idx = 0; idx < number_of_cells; ++idx)
            map[idx] = new MapCell{};

        // Read cells
        for (int idx = 0; idx < number_of_cells; ++idx) {
            auto curCell = map[idx];

            int type, neighbour;
            std::cin >> type >> curCell->resources;
            curCell->type = (MC_Type)type;
            curCell->idx = idx;

            for (int idx = 0; idx < 6; ++idx) {
                std::cin >> neighbour;

                if (neighbour != -1)
                    curCell->neighbours[idx] = map[neighbour];
            }
            std::cin.ignore();
        }

        // Read bases
        std::cin >> num_of_bases; std::cin.ignore();
        my_bases.resize(num_of_bases);
        opp_bases.resize(num_of_bases);

        for (int i = 0; i < num_of_bases; i++) {
            int my_base_index;
            std::cin >> my_base_index; std::cin.ignore();

            my_bases[i] = map[my_base_index];
        }
        for (int i = 0; i < num_of_bases; i++) {
            int opp_base_index;
            std::cin >> opp_base_index; std::cin.ignore();

            opp_bases[i] = map[opp_base_index];
        }
    }

    void ReadInput()
    {
        static bool isFirstCall = true;

        // Reset values
        my_ants = 0;
        opp_ants = 0;
        cur_crystals = 0;
        cur_eggs = 0;
        crystals.clear();
        eggs.clear();

        for (int i = 0; i < number_of_cells; i++) {
            int resources; // the current amount of eggs/crystals on this cell
            int my_ants_in_cell;   // the amount of your ants on this cell
            int opp_ants_in_cell;  // the amount of opponent ants on this cell
            std::cin >> resources >> my_ants_in_cell >> opp_ants_in_cell; std::cin.ignore();
            map[i]->resources = resources;
            map[i]->my_ants = my_ants_in_cell;
            map[i]->opp_ants = opp_ants_in_cell;
            my_ants += my_ants_in_cell;
            opp_ants += opp_ants_in_cell;


            if (isFirstCall)
            {
                if (map[i]->type == MC_CRYSTAL)
                    init_crystals += resources;
                else if (map[i]->type == MC_EGGS)
                    init_eggs += resources;
            }

            if (map[i]->type == MC_CRYSTAL && resources != 0) {
                cur_crystals += resources;
                crystals.emplace_back(map[i]);
            }
            else if (map[i]->type == MC_EGGS && resources != 0) {
                cur_eggs += resources;
                eggs.emplace_back(map[i]);
            }
        }

        isFirstCall = false;
    }

} GInf;

std::vector<MapCell*>& MapCell::map = GInf.map;

class Command {
    std::stringstream ss;
    bool bEmpty = true;

public:
    bool empty() const { return bEmpty; }

    void Clear() { ss.str(""); bEmpty = true; }
    void Submit()
    {
        std::cout << ss.str() << std::endl;
    }

    void AddBeacon(int idx, int strength)
    {
        bEmpty = false;
        ss << "BEACON " << idx << ' ' << strength << ';';
    }

    void AddLine(int idxSrc, int idxDest,int strength)
    {
        bEmpty = false;
        ss << "LINE " << idxSrc << ' ' << idxDest << ' ' << strength << ';';
    }

    void AddWait()
    {
        bEmpty = false;
        ss << "WAIT";
    }
} GCommand;

class CachebleWaveAlg
{
public:

    struct SearchResult {
        int distance;
        std::deque<int> path;
    };


    CachebleWaveAlg(CachebleWaveAlg const&) = delete;
    void operator=(CachebleWaveAlg const&) = delete;
    static CachebleWaveAlg& getInstance()
    {
        static CachebleWaveAlg instance;
        return instance;
    }

    SearchResult* waveAlgorithm(int srcIdx, int destIdx, bool bEggsEmptyOnly)
    {
        //DBG_MSG_V2("waveAlgorithm", srcIdx, destIdx);
        auto it = cache.find({ srcIdx, destIdx });
        if (it != cache.end())
            return &it->second;

        //it = cache.find({ destIdx, srcIdx });
        //if (it != cache.end())
        //{
        //    cache[{destIdx, srcIdx}] = it->second;
        //    return &it->second;
        //}

        int distanceToDest = std::numeric_limits<int>::max();
        std::vector<int> distToCell;
        std::deque<int> path, bestPath;
        distToCell.assign(GInf.number_of_cells, std::numeric_limits<int>::max());

        waveAlgorithm_impl(srcIdx, destIdx, distToCell, 0, distanceToDest, path, bestPath, bEggsEmptyOnly);
        if (bEggsEmptyOnly)
        {
            cache_eggs_empty_only[{srcIdx, destIdx}] = { distToCell[destIdx], std::move(bestPath) };
            return &cache_eggs_empty_only[{srcIdx, destIdx}];
        }
        else
        {
            cache[{srcIdx, destIdx}] = { distToCell[destIdx], std::move(bestPath) };
            return &cache[{srcIdx, destIdx}];
        }
    }

private:


    CachebleWaveAlg() {}

    void waveAlgorithm_impl(int curIdx, int destIdx, std::vector<int>& distToCell, int distance, int& distanceToDest, std::deque<int>& path, std::deque<int>& bestPath, bool bEggsEmptyOnly)
    {
        auto& neighbours = GInf.map[curIdx]->neighbours;

        if (distToCell[curIdx] <= distance || distanceToDest <= distance)
            return;

        if (bEggsEmptyOnly && GInf.map[curIdx]->type == MC_CRYSTAL && GInf.map[curIdx]->resources != 0)
            return;

        path.emplace_back(curIdx);
        distToCell[curIdx] = distance;

        // destIdx was found
        if (curIdx == destIdx) {
            if (distance < distanceToDest) {
                distanceToDest = distance;
                bestPath = path;
            }

            path.pop_back();
            return;
        }

        for (auto pNeighbour : neighbours)
            if (pNeighbour != nullptr)
                waveAlgorithm_impl(pNeighbour->idx, destIdx, distToCell, distance + 1, distanceToDest, path, bestPath, bEggsEmptyOnly);

        path.pop_back();
    }

    std::map<std::pair<int, int>, SearchResult> cache;
    std::map<std::pair<int, int>, SearchResult> cache_eggs_empty_only;
};

struct EggPath {
    MapCell* targetCell;
    int distance;
    int resource;
    std::deque<int> path;
};

struct IdxStrength {
    int idx;
    int strength;
};
using EvaluateSteps = std::vector<std::vector<IdxStrength>>;

class StateMachine {
public:
    StateMachine(std::vector<EggPath>& eggs, int antsOnBase) : m_eggs(eggs), m_totalAnts(antsOnBase)
    {}

    void EvaluateStep()
    {
        GreedyEvaluateAntsStep();
        MoveAnts();
        GainCrystals();
    }

    void GreedyEvaluateAntsStep()
    {
        auto it = std::find_if(m_eggs.begin(), m_eggs.end(), [](const auto& o) {return o.resource != 0; });
        if (it == m_eggs.end())
        {
            m_bStop = true;
            return;
        }

        int antsPerCell = m_totalAnts / (it->distance + 1);
        m_evalSteps.emplace_back(std::vector<IdxStrength>{});
        for (int j = 0; j < it->path.size(); ++j)
            m_evalSteps.back().emplace_back(IdxStrength{ it->path[j], antsPerCell});
    }

    void MoveAnts()
    {

    }

    void GainCrystals()
    {

    }

private:
    bool m_bStop = false;

    int m_totalAnts = 0;
    int m_minAntsOnPath = 0;
    std::vector<EggPath>& m_eggs;

    EvaluateSteps m_evalSteps;
};

struct CellNeighbours {
    MapCell* pCell = nullptr;
    std::vector<MapCell*> neighbours;
    int totalResources = 0;
    int distFromBase = 0;

    int getCurCellResources() { return std::accumulate(neighbours.begin(), neighbours.end(), pCell->resources, [](int u, auto pCell) { return u + pCell->resources; }); }
};

CachebleWaveAlg::SearchResult* getNearestCellToPos(std::vector<MapCell*>& cells, int posIdx, bool bEggsEmptyOnly)
{
    if (cells.empty())
        return nullptr;

    auto& cAlg = CachebleWaveAlg::getInstance();

    CachebleWaveAlg::SearchResult* pSr = nullptr;
    int nearestDist = std::numeric_limits<int>::max();
    for (auto pCell : cells)
    {
        if (pCell->resources == 0)
            continue;

        auto pSearchRes = cAlg.waveAlgorithm(posIdx, pCell->idx, bEggsEmptyOnly);
        if (pSearchRes->distance < nearestDist)
        {
            nearestDist = pSearchRes->distance;
            pSr = pSearchRes;
        }
    }

    return pSr;
}

CachebleWaveAlg::SearchResult* getFarAwayCellToPos(std::vector<MapCell*>& cells, int posIdx, bool bEggsEmptyOnly)
{
    if (cells.empty())
        return nullptr;

    auto& cAlg = CachebleWaveAlg::getInstance();

    CachebleWaveAlg::SearchResult* pSr = nullptr;
    int maxDist = 0;
    for (auto pCell : cells)
    {
        if (pCell->resources == 0)
            continue;

        auto pSearchRes = cAlg.waveAlgorithm(posIdx, pCell->idx, bEggsEmptyOnly);
        if (pSearchRes->distance > maxDist)
        {
            maxDist = pSearchRes->distance;
            pSr = pSearchRes;
        }
    }

    return pSr;
}

CellNeighbours GetBestCellSpot(std::vector<MapCell*> cells, int maxSpotDistance)
{
    auto& cAlg = CachebleWaveAlg::getInstance();
    CellNeighbours ret;

    // Find all spots
    std::vector<CellNeighbours> neighbourCrystalsMap;
    for (auto pCrystal1 : cells)
    {
        auto pSearchRes = cAlg.waveAlgorithm(GInf.my_bases[0]->idx, pCrystal1->idx, false);

        neighbourCrystalsMap.emplace_back(CellNeighbours{ pCrystal1, {}, pCrystal1->resources, pSearchRes->distance + 1 });
        for (auto pCrystal2 : cells)
        {
            if (pCrystal1 == pCrystal2)
                continue;

            auto pSearchRes = cAlg.waveAlgorithm(pCrystal1->idx, pCrystal2->idx, false);
            if (pSearchRes->distance < maxSpotDistance) {
                neighbourCrystalsMap.back().neighbours.emplace_back(pCrystal2);
                neighbourCrystalsMap.back().totalResources += pCrystal2->resources;
            }
        }

        if (neighbourCrystalsMap.back().neighbours.empty())
            neighbourCrystalsMap.pop_back();
    }

    if (!neighbourCrystalsMap.empty())
    {
        // Select best spot
        std::sort(neighbourCrystalsMap.begin(), neighbourCrystalsMap.end(), [](const auto& lhd, const auto& rhd) {
            return lhd.totalResources / lhd.distFromBase > rhd.totalResources / rhd.distFromBase;
            });

        ret = neighbourCrystalsMap[0];
    }

    return ret;
}

class Strategy
{
public:
    void Prepare()
    {
        if (m_curCrystalSpot.pCell != nullptr && m_curCrystalSpot.getCurCellResources() == 0)
            m_curCrystalSpot.pCell = nullptr;

        if (m_curEggSpot.pCell != nullptr && m_curEggSpot.getCurCellResources() == 0)
            m_curEggSpot.pCell = nullptr;

        if (m_curBestEnemyCrystal != -1 && GInf.map[m_curBestEnemyCrystal]->resources == 0)
        {
            m_curBestEnemyCrystal = -1;
        }

        if (m_curBestCrystal != -1 && GInf.map[m_curBestCrystal]->resources == 0)
        {
            m_curBestCrystal = -1;
        }
    }

    void Execute()
    {
        DBG_MSG_V2("Strategy::Execute+", m_curBestEnemyCrystal, m_curBestCrystal);


        auto& cAlg = CachebleWaveAlg::getInstance();
        if (m_curStep == 0)
        {
            auto pSearchRes = cAlg.waveAlgorithm(GInf.my_bases[0]->idx, GInf.opp_bases[0]->idx, false);
            m_mapDiameter = pSearchRes->distance + 1;

            for (auto pCell1 : GInf.map)
                for (auto pCell2 : GInf.map)
                {
                    auto pSearchRes = cAlg.waveAlgorithm(pCell1->idx, pCell2->idx, false);
                    m_mapDiameter = std::max(pSearchRes->distance + 1, m_mapDiameter);
                }
        }


        std::vector<EggPath> eggsSR;
        for (auto pEgg : GInf.eggs)
        {
            if (pEgg->resources == 0)
                continue;

            auto pSearchRes = cAlg.waveAlgorithm(GInf.my_bases[0]->idx, pEgg->idx, false);
            //DBG_MSG_V2("X0", pSearchRes->distance, pSearchRes->path.size());
            if (pSearchRes) eggsSR.emplace_back(EggPath{ pEgg, pSearchRes->distance, pEgg->resources, pSearchRes->path });
            else            DBG_MSG_S("ERR: pSearchRes == nullptr");
        }

        std::sort(eggsSR.begin(), eggsSR.end(), [](const auto& lhd, const auto& rhd) { return lhd.distance < rhd.distance; });
        eggsSR.erase(
            std::remove_if(eggsSR.begin(), eggsSR.end(), [](const auto& o) { return o.distance == std::numeric_limits<int>::max() || o.resource == 0; }),
            eggsSR.end());

        if (eggsSR.empty() || (!eggsSR.empty() && m_mapDiameter / 2.2 < eggsSR[0].distance))
        {
            DBG_MSG_S("crystal");

            // Find best crystals spot
            if (m_curCrystalSpot.pCell == nullptr)
            {
                static const int MAX_SPOT_DISTANCE = 3;

                m_curCrystalSpot = GetBestCellSpot(GInf.crystals, GInf.MAX_CRYSTALS_SPOT_DISTANCE);
                /*
                std::vector<CellNeighbours> neighbourCrystalsMap;
                for (auto pCrystal1 : GInf.crystals)
                {
                    auto pSearchRes = cAlg.waveAlgorithm(GInf.my_bases[0]->idx, pCrystal1->idx, false);

                    neighbourCrystalsMap.emplace_back(CellNeighbours{ pCrystal1, {}, pCrystal1->resources, pSearchRes->distance + 1 });
                    for (auto pCrystal2 : GInf.crystals)
                    {
                        if (pCrystal1 == pCrystal2)
                            continue;

                        auto pSearchRes = cAlg.waveAlgorithm(pCrystal1->idx, pCrystal2->idx, false);
                        if (pSearchRes->distance < MAX_SPOT_DISTANCE) {
                            neighbourCrystalsMap.back().neighbours.emplace_back(pCrystal2);
                            neighbourCrystalsMap.back().totalResources += pCrystal2->resources;
                        }
                    }

                    if (neighbourCrystalsMap.back().neighbours.empty())
                        neighbourCrystalsMap.pop_back();
                }

                if (!neighbourCrystalsMap.empty())
                {
                    // Select best spot
                    std::sort(neighbourCrystalsMap.begin(), neighbourCrystalsMap.end(), [](const auto& lhd, const auto& rhd) {
                        return lhd.totalResources / lhd.distFromBase > rhd.totalResources / rhd.distFromBase;
                        });

                    m_curCrystalSpot = neighbourCrystalsMap[0];
                }
                */
            }

            // If crystal spot exists - prefer it, otherwise - just get nearest crystal
            if (m_curCrystalSpot.pCell != nullptr)
            {
                auto& cAlg = CachebleWaveAlg::getInstance();
                auto& n = m_curCrystalSpot.neighbours;

                n.erase(std::remove_if(n.begin(), n.end(), [](auto pCell) { return pCell->resources == 0; }), n.end());

                DBG_MSG_S("=============== Evaluate crystal spot logic ===============");
                DBG_V2(m_curCrystalSpot.pCell->idx, m_curCrystalSpot.neighbours.size());


                int totalDistance = std::accumulate(n.begin(), n.end(), m_curCrystalSpot.distFromBase, [this, &cAlg](int u, MapCell* pCell) {
                    auto pSearchRes = cAlg.waveAlgorithm(m_curCrystalSpot.pCell->idx, pCell->idx, false);

                    return u + pSearchRes->distance;
                    });

                int antsPerCell = GInf.my_ants / totalDistance;
                GCommand.AddLine(GInf.my_bases[0]->idx, m_curCrystalSpot.pCell->idx, antsPerCell);

                for (auto pNeighbour : m_curCrystalSpot.neighbours)
                {
                    auto pSearchRes = cAlg.waveAlgorithm(m_curCrystalSpot.pCell->idx, pNeighbour->idx, false);
                    GCommand.AddLine(pSearchRes->path[1], pNeighbour->idx, antsPerCell);
                }
            }
            else
            {
                std::vector<CachebleWaveAlg::SearchResult*> crystalsSR;
                for (auto pCrystal : GInf.crystals)
                {
                    if (pCrystal->resources == 0)
                        continue;

                    auto pSearchRes = cAlg.waveAlgorithm(GInf.my_bases[0]->idx, pCrystal->idx, false);
                    DBG_MSG_V3("X0 crystals", pCrystal->idx, pSearchRes->distance, pSearchRes->path.size());
                    if (pSearchRes) crystalsSR.emplace_back(pSearchRes);
                    else            DBG_MSG_S("ERR: pSearchRes == nullptr");
                }

                auto enemyCrystalsSR = crystalsSR; // create a copy

                //// Best enemy crystal
                //if (m_curBestEnemyCrystal == -1)
                //{
                //    DBG_MSG_S("=============== Find new 'best enemy crystal' logic ===============");
                //    enemyCrystalsSR.erase(
                //        std::remove_if(enemyCrystalsSR.begin(), enemyCrystalsSR.end(), [this](const auto pO) { return pO->distance < m_mapDiameter / 2; }),
                //        enemyCrystalsSR.end());
                //
                //
                //    //DBG_MSG_S("enemyCrystalsSR:");
                //    //for (auto pCell : enemyCrystalsSR)
                //    //    std::cerr << "back: " << pCell->path.back() << ' ' << "size: " << pCell->path.size() << " first " << pCell->path.front();
                //    //std::cerr << '\n';
                //
                //    std::sort(enemyCrystalsSR.begin(), enemyCrystalsSR.end(), [](const auto lhd, const auto rhd) { return GInf.map[lhd->path.back()]->resources / (lhd->distance + 1) > GInf.map[rhd->path.back()]->resources / (rhd->distance + 1); });
                //
                //
                //    if (!enemyCrystalsSR.empty())
                //    {
                //        auto pNearestCrystalSR = enemyCrystalsSR[0];
                //        m_curBestEnemyCrystal = pNearestCrystalSR->path.back();
                //
                //        if (m_curBestEnemyCrystal == GInf.my_bases[0]->idx)
                //            m_curBestEnemyCrystal = pNearestCrystalSR->path.front();
                //    }
                //}
                //if (m_curBestEnemyCrystal != -1)
                //{
                //    DBG_MSG_S("=============== Evaluate 'best enemy crystal' logic ===============");
                //    DBG_V2(GInf.my_bases[0]->idx, m_curBestEnemyCrystal);
                //    GCommand.AddLine(GInf.my_bases[0]->idx, m_curBestEnemyCrystal, 1);
                //}

                // Best crystal
                if (m_curBestCrystal == -1)
                {
                    DBG_MSG_S("=============== Find new 'best crystal' logic ===============");
                    std::sort(crystalsSR.begin(), crystalsSR.end(), [](const auto lhd, const auto rhd) { return GInf.map[lhd->path.back()]->resources / (lhd->distance + 1) > GInf.map[rhd->path.back()]->resources / (rhd->distance + 1); });

                    if (!crystalsSR.empty())
                    {
                        int antsPerCell = GInf.my_ants;
                        if (!enemyCrystalsSR.empty())
                            antsPerCell -= (enemyCrystalsSR[0]->distance + 1)*2;

                        antsPerCell /= (crystalsSR[0]->distance + 1);
                        GCommand.AddLine(GInf.my_bases[0]->idx, crystalsSR[0]->path.back(), antsPerCell);
                        m_curBestCrystal = crystalsSR[0]->path.back();
                        m_currFarAwayCrystalStrengts = antsPerCell;
                    }
                    else
                    {
                        DBG_MSG_S("ERROR: crystalsSR.empty() == true");
                    }
                }
                if (m_curBestCrystal != -1)
                {
                    DBG_MSG_S("=============== Evaluate 'best crystal' logic ===============");
                    DBG_V2(GInf.my_bases[0]->idx, m_curBestCrystal);
                    GCommand.AddLine(GInf.my_bases[0]->idx, m_curBestCrystal, m_currFarAwayCrystalStrengts);
                }

                //for (auto pCell : GInf.crystals)
                //{
                //    if (pCell->resources == 0 || pCell->idx == m_curBestCrystal || pCell->idx == m_curBestEnemyCrystal)
                //        continue;
                //
                //    GCommand.AddLine(GInf.my_bases[0]->idx, pCell->idx, 1);
                //}
            }
        }
        else
        {
            DBG_MSG_S("eggs");

            if (m_curEggSpot.pCell == nullptr)
            {
                auto eggsNearestToMe = GInf.eggs;
                eggsNearestToMe.erase(std::remove_if(eggsNearestToMe.begin(), eggsNearestToMe.end(), [&cAlg](auto pCell) {
                    auto pSearchResMy = cAlg.waveAlgorithm(GInf.my_bases[0]->idx, pCell->idx, false);
                    auto pSearchResOpp = cAlg.waveAlgorithm(GInf.opp_bases[0]->idx, pCell->idx, false);

                    return pSearchResOpp->distance < pSearchResMy->distance;
                    }),
                    eggsNearestToMe.end());

                if (!eggsNearestToMe.empty())
                {
                    m_curEggSpot = GetBestCellSpot(eggsNearestToMe, GInf.MAX_EGGS_SPOT_DISTANCE);
                }
            }

            if (m_curEggSpot.pCell != nullptr)
            {
                auto& cAlg = CachebleWaveAlg::getInstance();
                auto& n = m_curEggSpot.neighbours;

                DBG_MSG_S("=============== Evaluate eggs spot logic ===============");
                DBG_V2(m_curEggSpot.pCell->idx, m_curEggSpot.neighbours.size());

                int totalDistance = std::accumulate(n.begin(), n.end(), m_curEggSpot.distFromBase, [this, &cAlg](int u, MapCell* pCell) {
                    auto pSearchRes = cAlg.waveAlgorithm(m_curEggSpot.pCell->idx, pCell->idx, false);

                    return u + pSearchRes->distance;
                    });

                int antsPerCell = GInf.my_ants / totalDistance;
                GCommand.AddLine(GInf.my_bases[0]->idx, m_curEggSpot.pCell->idx, antsPerCell);

                for (auto pNeighbour : m_curEggSpot.neighbours)
                {
                    auto pSearchRes = cAlg.waveAlgorithm(m_curEggSpot.pCell->idx, pNeighbour->idx, false);
                    GCommand.AddLine(pSearchRes->path[1], pNeighbour->idx, antsPerCell);
                }
            }
            else
            {
                DBG_MSG_S("=============== Evaluate simple eggs logic ===============");
                DBG_V3(GInf.my_ants, eggsSR[0].targetCell->idx, eggsSR[0].distance);

                int antsPerCell = GInf.my_ants / (eggsSR[0].distance + 1);
                auto eggCell = eggsSR[0].targetCell;
                if (antsPerCell <= eggCell->resources)
                {
                    DBG_V(antsPerCell);
                    GCommand.AddLine(eggsSR[0].path.front(), eggsSR[0].path.back(), antsPerCell);
                }
                else
                {
                    int antsPerCell_0 = eggCell->resources;
                    int antsPerCell_1 = (GInf.my_ants - eggCell->resources * eggsSR[0].distance) / (eggsSR[1].distance + 1);
                    DBG_V3(GInf.my_ants, antsPerCell_0, antsPerCell_1);

                    if (eggsSR.size() == 1)
                        return; // no more eggs

                    //for (int idx : eggsSR[1].path)
                    //    GCommand.AddBeacon(idx, antsPerCell_1);
                    //
                    //for (int idx : eggsSR[0].path)
                    //    GCommand.AddBeacon(idx, antsPerCell_0);
                    GCommand.AddLine(eggsSR[1].path.front(), eggsSR[1].path.back(), antsPerCell);
                    GCommand.AddLine(eggsSR[0].path.front(), eggsSR[0].path.back(), antsPerCell);
                }
            }
        }


        DBG_MSG_V2("Strategy::Execute-", m_curBestEnemyCrystal, m_curBestCrystal);
        m_curStep++;
    }

private:
    int m_mapDiameter = 0;
    int m_curStep = 0;
    int m_evalStepsLen = 0;
    //EggsCollectorAlg::EvaluateSteps m_evalSteps;

    int m_curBestEnemyCrystal = -1;
    int m_currFarAwayCrystalStrengts = -1;
    int m_curBestCrystal = -1;

    CellNeighbours m_curCrystalSpot;
    CellNeighbours m_curEggSpot;
};

int main()
{
    Strategy strategy{};

    // game loop
    GInf.Initialize();
    while (1) {
        GCommand.Clear();
        GInf.ReadInput();

        strategy.Prepare();
        strategy.Execute();

        if (GCommand.empty())
            GCommand.AddWait();

        GCommand.Submit();
    }
}