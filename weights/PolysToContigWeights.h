#ifndef __GEODA_CENTER_POLYS_TO_CONTIG_WEIGHTS_H__
#define __GEODA_CENTER_POLYS_TO_CONTIG_WEIGHTS_H__

#include <cfloat>

#include "../geofeature.h"
#include "GalWeight.h"

#define CMP_DBL_EPSILON sqrt(DBL_EPSILON)

static const int GdaConst_EMPTY = -1;

/**
 BasePartition
 */
class BasePartition  {
    protected :
    int         elements, cells;
    int *       cell;
    int *       next;
    double      step;
    public :
    BasePartition(const int els= 0, const int cls= 0, const double range= 0);
    virtual ~BasePartition();
    void virtual alloc(const int els, const int cls, const double range);
    int         Cells() const  {  return cells;  };
    double      Step() const  {  return step;  };
    virtual void include(const int incl, const double range)  {
        int where = (int) floor(range/step);
        // if (where < -1 || where > cells || incl < 0 || incl >= elements)
        //     cout << " BasePartition: incl= " << incl << " location= "
        //          << where << " els= " << elements << " cells= "
        //          << cells << endl;
        if (where < 0) where= 0;
        else if (where >= cells) where= cells-1;
        next [ incl ] = cell [ where ];
        cell [ where ] = incl;
        return;
    };
    
    int first(const int cl) const  {  return cell [ cl ];  };
    int tail(const int elt) const  {  return next [ elt ];  };
};


/**
 PartitionP
 */
class PartitionP : public BasePartition  {
    private :
    int *       cellIndex;
    int *       previous;
    public :
    PartitionP(const int els= 0, const int cls= 0, const double range= 0);
    ~PartitionP();
    void alloc(const int els, const int cls, const double range);

    using BasePartition::include; // tell the compiler we want both the include from Base and ours
    virtual void include(const int incl);
    void initIx(const int incl, const double range)  {
        int cl= (int) floor(range / step);
        // if (cl < -1 || cl > cells || incl < 0 || incl >= elements)
        //     cout << "PartitionP: incl= " << incl << " at " << cl << endl;
        if (cl < 0) cl= 0;
        else if (cl >= cells) cl= cells-1;
        cellIndex[ incl ] = cl;
        return;
    };
    int inTheRange(const double range) const
    {
        if (range < 0 || range/step > cells + CMP_DBL_EPSILON) return -1;
        int where= (int) floor(range / step);
        if (where < 0) where= 0;
        else if (where >= cells) --where;
        return where;
    }
    void remove(const int del);
    void cleanup(const BasePartition &p, const int cl)  {
        for (int cnt= p.first(cl); cnt != GdaConst_EMPTY; cnt= p.tail(cnt))
            remove(cnt);
    }
};

/*
 * Polygon Partition
 */
class PolygonPartition
{
    protected :
    gda::PolygonContents     *poly;
    
    BasePartition       pX;
    PartitionP          pY;
    int *               nbrPoints;
    
    int prev(const int pt) const
    {
        int ix= nbrPoints[pt];
        return (ix >= 0) ? pt-1 : -ix;
    }
    int succ(const int pt) const
    {
        int ix= nbrPoints[pt];
        return (ix >= 0) ? ix : pt+1;
    }
    
    public :
    int                 NumPoints;
    int                 NumParts;
    
    PolygonPartition(gda::PolygonContents* _poly)
    : pX(), pY(), nbrPoints(NULL) {
        poly = _poly;
        NumPoints = poly->num_points;
        NumParts = poly->num_parts;
    }
    ~PolygonPartition();
    
    gda::Point* GetPoint(const int i){ return &poly->points[i];}
    int GetPart(int i){ return (int)poly->parts[i]; }
    double GetMinX(){ return (double)poly->box[0]; }
    double GetMinY(){ return (double)poly->box[1]; }
    double GetMaxX(){ return (double)poly->box[2]; }
    double GetMaxY(){ return (double)poly->box[3]; }
    
    int  MakePartition(int mX= 0, int mY= 0);
    void MakeSmallPartition(const int mX, const double Start,
                            const double Stop);
    void MakeNeighbors();
    bool edge(PolygonPartition &p, const int host, const int guest,
              double precision_threshold);
    int sweep(PolygonPartition & guest, bool is_queen,
              double precision_threshold=0.0);
};


typedef struct Ref  {
    int next, prev;
    Ref(const int nxt= -1, const int prv= -1) : next(nxt), prev(prv) {};
} RefStruct;

typedef RefStruct* RefPtr;

/*
 PartitionM
 */
class PartitionM  {
    private :
    double      step;
    int         elements, cells;
    int *       cell;
    int *       cellIndex;
    int *       lastIndex;
    RefPtr *    Refs;
    public :
    PartitionM(const int els, const int cls, const double range);
    virtual ~PartitionM();
    
    void include(const int incl);
    void remove(const int del);
    void initIx(const int incl, const double lwr, const double upr);
    int lowest (const int el) const  {  return cellIndex [ el ];  };
    int upmost(const int el) const  {  return lastIndex [ el ];  };
    int first(const int cl) const  {  return cell[ cl ];  };
    int tail(const int elt, const int cl) const  {
        return Refs[elt][cl-cellIndex[elt]].next;
    };
    int Sum() const;
};

GalElement* PolysToContigWeights(gda::MainMap& main,
                                 bool is_queen,
                                 double precision_threshold=0.0);

#endif
