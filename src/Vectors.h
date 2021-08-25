#ifndef VECTORS_H
#define VECTORS_H
/* Create a header for defining vectors and vector operations */
#include <cmath>
#include <limits>

template<typename T, unsigned DIM> 
class vec
{
    static_assert(DIM > 0, "Dimension must be greater than 0");

    public:
        template<typename... Args>
        vec(Args&&... args) : data{T(args)...}
        {
            if(sizeof...(Args) != 1)
                static_assert(sizeof...(Args) == DIM, "Invalid number of constructor arguments.");
        }

    private:
        T data[DIM];

};

/*********************************************************/
/*****************  2D VECTOR CLASS **********************/
/*********************************************************/
template <typename T>
class vec<T,2>
{
    public:
    
        vec<T,2>(){}

        vec<T,2>(T const a): data{a,a} {}

        vec<T,2>(T const a, T const b): data{a,b} {}

        vec<T,2>(vec<T,2> const& a) {*this = a;}
        
        void zero()
        {
            data[0] = 0.0; data[1] = 0.0;
        }

        static vec<T,2> Zero() 
        {
            return vec<T,2>(0.0,0.0);
        }

        inline T dot(vec<T,2> const& a) const
        {
            return data[0]*a[0] + data[1]*a[1];
        }

        inline T cross(vec<T,2> const& a) const
        {
            return (data[0]*a[1] - data[1]*a[0]);
        }
        
        inline T const norm() const
        {
            return sqrt(data[0]*data[0] + data[1]*data[1]);
        }

        inline T squaredNorm() const
        {
            return (data[0]*data[0] + data[1]*data[1]);
        }

        inline vec<T,2> normalized() const
        {
            return vec<T,2>(data[0]/ norm(), data[1] / norm());
        }

        /* Standard operators */
        inline vec<T,2> operator + (vec<T,2> const& a) const
        {
            return vec(data[0]+a[0], data[1]+a[1]);
        }

        inline vec<T,2> operator + (T const a) const
        {
            return vec<T,2>(data[0] + a, data[1] + a);
        }

        inline vec<T,2> operator - (vec<T,2> const& a) const
        {
            return vec(data[0]-a[0], data[1]-a[1]);
        }

        inline vec<T,2> operator - (T const a) const
        {
            return vec<T,2>(data[0] - a, data[1] - a);
        }

        inline vec<T,2> operator * (T const a) const
        {
            return vec<T,2>(data[0] * a, data[1] * a);
        }

        inline vec<T,2> operator / (T const a) const
        {
            return vec<T,2>(data[0] / a, data[1] / a);
        }

        /* = operators */
        vec<T,2>& operator += (vec<T,2> const& a)
        {
            data[0] = data[0] + a[0];
            data[1] = data[1] + a[1];
            return *this;
        }

        vec<T,2>& operator += (T const a)
        {
            data[0] = data[0] + a;
            data[1] = data[1] + a;
            return *this;
        }

        vec<T,2>& operator -= (vec<T,2> const& a)
        {
            data[0] = data[0] - a[0];
            data[1] = data[1] - a[1];
            return *this;
        }

        vec<T,2>& operator -= (T const a)
        {
            data[0] = data[0] - a;
            data[1] = data[1] - a;
            return *this;
        }

        vec<T,2>& operator *= (T const a)
        {
            data[0] = data[0] * a;
            data[1] = data[1] * a;
            return *this;
        }

        vec<T,2>& operator /= (T const a)
        {
            data[0] = data[0] / a;
            data[1] = data[1] / a;
            return *this;
        }

        vec<T,2>& operator = (vec<T,2> const& a)
        {
            data[0] = a[0];
            data[1] = a[1];
            return *this;
        }

        /* Equality operators */
        bool operator == (vec<T,2> const& a) const
        {
            return (data[0] == a[0] && data[1] == a[1]);
        }

        bool operator != (vec<T,2> const& a) const
        {
            return (data[0] != a[0] || data[1] != a[1]);
        }

        // bool operator >= (vec<T,2> const& a) const
        // {
        //     return (data[0] >= a[0] && data[1] >= a[1]);
        // }

        // bool operator <= (vec<T,2> const& a) const
        // {
        //     return (data[0] <= a[0] && data[1] <= a[1]);
        // }

        // bool operator > (vec<T,2> const& a) const
        // {
        //     return (data[0] > a[0] && data[1] > a[1]);
        // }

        // bool operator < (vec<T,2> const& a) const
        // {
        //     return (data[0] < a[0] && data[1] < a[1]);
        // }
        
        T& operator [](size_t const index) 
        {
            return data[index];
        }

        T const& operator [](size_t const index) const
        {
            return data[index];
        }

        T& operator ()(size_t const index)
        {
            return data[index];
        }

        T const& operator ()(size_t const index) const
        {
            return data[index];
        }

        size_t const size() const {return 2;}

    private:
        T data[2];
};

template<typename T>
inline vec<T,2> const operator*(T const& a, vec<T,2> const& b)
{
    return b*a;
}

template<typename T>
inline vec<T,2> const operator-(vec<T,2> const& b)
{
    return -b;
}

/*********************************************************/
/*****************  3D VECTOR CLASS **********************/
/*********************************************************/
template <typename T>
class vec<T,3>
{
    public:
        vec<T,3>(){}

        vec<T,3>(T const a): data{a,a,a} {}

        vec<T,3>(T const a, T const b, T const c): data{a,b,c} {}

        vec<T,3>(vec<T,3> const& a) {*this = a;}

        void zero()
        {
            data[0] = 0.0; data[1] = 0.0; data[2] = 0.0;
        }

        static vec<T,3> Zero() 
        {
            return vec<T,3>(0.0,0.0,0.0);
        }

        inline T dot(vec<T,3> const& a) const
        {
            return data[0]*a[0] + data[1]*a[1] + data[2]*a[2];
        }

        inline vec<T,3> cross(vec<T,3> const& a) const
        {
            return 
            vec<T,3>((data[1]*a[2] - data[2]*a[1]),
                (data[2]*a[0] - data[0]*a[2]),
                (data[0]*a[1] - data[1]*a[0]));
        }

        inline T const norm() const
        {
            return sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
        }

        inline T const squaredNorm() const
        {
            return (data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
        }

        inline vec<T,3> normalized() const
        {
            return *this / norm();
        }

        /* Standard operators */
        inline vec<T,3> operator+(vec<T,3> const a) const
        {
            return vec(data[0] + a[0], data[1] + a[1], data[2] + a[2]);
        }

        inline vec<T,3> operator+(T const a) const
        {
            return vec(data[0] + a, data[1] + a, data[2] + a);
        }

        inline vec<T,3> operator - (vec<T,3> const a) const
        {
            return vec<T,3>(data[0] - a[0], data[1] - a[1], data[2] - a[2]);
        }

        inline vec<T,3> operator - (T const a) const
        {
            return vec<T,3>(data[0] - a, data[1] - a, data[2] - a);
        }

        inline vec<T,3> operator / (T const a) const
        {
            return vec<T,3>(data[0] / a, data[1] / a, data[2] / a);
        }

        inline vec<T,3> operator * (T const a) const
        {
            return vec<T,3>(data[0] * a, data[1] * a, data[2] * a);
        }

        vec<T,3>& operator +=(vec<T,3> const& a)
        {
            data[0] = data[0] + a[0];
            data[1] = data[1] + a[1];
            data[2] = data[2] + a[2];
            return *this;
        }

        /* = operators */
        vec<T,3>& operator += (T const& a)
        {
            data[0] = data[0] + a;
            data[1] = data[1] + a;
            data[2] = data[2] + a;
            return *this;
        }

        vec<T,3>& operator -= (vec<T,3> const& a)
        {
            data[0] = data[0] - a[0];
            data[1] = data[1] - a[1];
            data[2] = data[2] - a[2];
            return *this;
        }

        vec<T,3>& operator -=(T const& a)
        {
            data[0] = data[0] - a;
            data[1] = data[1] - a;
            data[2] = data[2] - a;
            return *this;
        }

        vec<T,3>& operator /= (T const a)
        {
            data[0] = data[0] / a;
            data[1] = data[1] / a;
            data[2] = data[2] / a;
            return *this;
        }

        vec<T,3>& operator *= (T const a)
        {
            data[0] = data[0] * a;
            data[1] = data[1] * a;
            data[2] = data[2] * a;
            return *this;
        }

        vec<T,3>& operator = (vec<T,3> const& a)
        {
            data[0] = a[0];
            data[1] = a[1];
            data[2] = a[2];
            return *this;
        }

        /* Equality operators */
        bool operator == (vec<T,3> const& a) const
        {
            return (data[0] == a[0] && data[1] == a[1]  && data[2] == a[2]);
        }

        bool operator != (vec<T,3> const& a) const
        {
            return (data[0] != a[0] || data[1] != a[1] || data[2] != a[2]);
        }

        T& operator [](size_t const index) 
        {
            return data[index];
        }

        T const& operator [](size_t const index) const
        {
            return data[index];
        }

        T& operator ()(size_t const index)
        {
            return data[index];
        }

        T const& operator ()(size_t const index) const
        {
            return data[index];
        }

        size_t const size() const {return 3;} 

    private:
        T data[3];
};

template<typename T>
vec<T,3> const operator*(T const& a, vec<T,3> const& b)
{
    return b*a;
}

template<typename T>
vec<T,3> const operator-(vec<T,3> const& b)
{
    return vec<T,3>(-b[0],-b[1],-b[2]);
}

/*********************************************************/
/*****************  4D VECTOR CLASS **********************/
/*********************************************************/
template <typename T>
class vec<T,4>
{
    public:
        vec<T,4>(){}

        vec<T,4>(T const a): data{a,a,a,a} {}

        vec<T,4>(T const a, T const b, T const c, T const d): data{a,b,c,d} {}

        vec<T,4>(vec<T,4> const& a) {*this = a;}

        vec<T,4>(vec<T,3> const& a, T const b)
        {
            data[0] = a[0]; data[1] = a[1]; data[2] = a[2]; data[3] = b;
        }

        void zero()
        {
            data[0] = 0.0; data[1] = 0.0; data[2] = 0.0; data[3] = 0.0;
        }

        static vec<T,4> Zero() 
        {
            return vec<T,4>(0.0,0.0,0.0,0.0);
        }

        inline T dot(vec<T,4> const& a)
        {
            return data[0]*a[0] + data[1]*a[1] + data[2]*a[2] + data[3]*a[3];
        }

        /* Doesn't Tly exist in 4D, and not useful */
        // vec4 pcross(vec4 const& a)
        // {
        //     return 
        //     vec4((data[1]*a[2] - data[2]*a[1]),
        //          (data[2]*a[0] - data[0]*a[2]),
        //          (data[0]*a[1] - data[1]*a[0]),
        //           data[0]*a[1] - data[1]*a[0]));
        // }

        inline T const norm() const
        {
            return sqrt(data[0]*data[0] + data[1]*data[1] + 
                    data[2]*data[2] + data[3]*data[3]);
        }

        inline T const squaredNorm() const
        {
            return (data[0]*data[0] + data[1]*data[1] + 
                    data[2]*data[2] + data[3]*data[3]);
        }

        inline vec<T,4> normalized() const
        {
            return *this / norm();
        }

        /* Standard operators */
        inline vec<T,4> operator+(vec<T,4> const a) const
        {
            return vec<T,4>(data[0] + a[0], data[1] + a[1], 
                            data[2] + a[2], data[3] + a[3]);
        }

        inline vec<T,4> operator+(T const a) const
        {
            return vec<T,4>(data[0] + a, data[1] + a, 
                            data[2] + a, data[3] + a);
        }

        inline vec<T,4> operator - (vec<T,4> const a) const
        {
            return vec<T,4>(data[0] - a[0], data[1] - a[1], 
                            data[2] - a[2], data[3] - a[3]);
        }

        inline vec<T,4> operator - (T const a) const
        {
            return vec<T,4>(data[0] - a, data[1] - a, 
                            data[2] - a, data[3] - a);
        }

        inline vec<T,4> operator / (T const a) const
        {
            return vec<T,4>(data[0] / a, data[1] / a, 
                            data[2] / a, data[3] / a);
        }

        inline vec<T,4> operator * (T const a) const
        {
            return vec<T,4>(data[0] * a, data[1] * a, 
                            data[2] * a, data[3] * a);
        }

        vec<T,4>& operator +=(vec<T,4> const& a)
        {
            data[0] = data[0] + a[0];
            data[1] = data[1] + a[1];
            data[2] = data[2] + a[2];
            data[3] = data[3] + a[3];
            return *this;
        }

        /* = operators */
        vec<T,4>& operator += (T const& a)
        {
            data[0] = data[0] + a;
            data[1] = data[1] + a;
            data[2] = data[2] + a;
            data[3] = data[3] + a;
            return *this;
        }

        vec<T,4>& operator -= (vec<T,4> const& a)
        {
            data[0] = data[0] - a[0];
            data[1] = data[1] - a[1];
            data[2] = data[2] - a[2];
            data[3] = data[3] - a[3];
            return *this;
        }

        vec<T,4>& operator -=(T const& a)
        {
            data[0] = data[0] - a;
            data[1] = data[1] - a;
            data[2] = data[2] - a;
            data[3] = data[3] - a;
            return *this;
        }

        vec<T,4>& operator /= (T const a)
        {
            data[0] = data[0] / a;
            data[1] = data[1] / a;
            data[2] = data[2] / a;
            data[3] = data[3] / a;
            return *this;
        }

        vec<T,4>& operator *= (T const a)
        {
            data[0] = data[0] * a;
            data[1] = data[1] * a;
            data[2] = data[2] * a;
            data[3] = data[3] * a;
            return *this;
        }

        vec<T,4>& operator = (vec<T,4> const& a)
        {
            data[0] = a[0];
            data[1] = a[1];
            data[2] = a[2];
            data[3] = a[3];
            return *this;
        }

        /* Equality operators */
        bool operator == (vec<T,4> const& a) const
        {
            return (data[0] == a[0] && data[1] == a[1] && data[2] == a[2] && data[3] == a[3]);
        }

        bool operator != (vec<T,4> const& a) const
        {
            return (data[0] != a[0] || data[1] != a[1] || data[2] != a[2] || data[3] != a[3]);
        }

        T& operator [](size_t const index)
        {
            return data[index];
        }

        T const& operator [](size_t const index) const
        {
            return data[index];
        }

        T& operator ()(size_t const index)
        {
            return data[index];
        }

        T const& operator ()(size_t const index) const
        {
            return data[index];
        }

        size_t const size() const {return 4;}

    private:
        T data[4];
};

template<typename T>
vec<T,4> const operator*(T const& a, vec<T,4> const& b)
{
    return b*a;
}

template<typename T>
vec<T,4> const operator-(vec<T,4> const& b)
{
    return vec<T,4>(-b[0],-b[1],-b[2],-b[3]);
}

/************************************************************/
/*****************  MATRIX DEFINITIONS **********************/
/************************************************************/

/* Only consider square matrices */
template<typename T, unsigned DIM> 
class matrix
{
    static_assert(DIM > 0, "Dimension must be greater than 0");

    public:

    matrix(void) {}
    
    template<typename... Args>
    matrix(Args&&... args) : data{T(args)...} {}

    T const& operator [](size_t const index)  const
    {
        return data[index];
    }

    private:
        T data[DIM*DIM];
};

template<typename T, unsigned DIM>
matrix<T,DIM> const operator*(T const& a, matrix<T,DIM> const& b)
{
    return b*a;
}

template<typename T, unsigned DIM>
matrix<T,DIM> const operator-(matrix<T,DIM> const& b)
{
    return -b;
}

/*********************************************************/
/*****************  2D MATRIX CLASS **********************/
/*********************************************************/
template <typename T>
class matrix<T,2>
{
    public:
        matrix<T,2>(void) {}

        matrix<T,2>(T const a): data{a,a,a,a} {}

        matrix<T,2>(T const a, T const b, T const c, T const d): data{a,b,c,d} {}

        matrix<T,2>(matrix<T,2> const& a) {*this = a;}
        
        T const determinant() const
        {
            return (data[0]*data[3] - data[1]*data[2]);
        }

        matrix<T,2> const transpose()
        {
            return matrix<T,2>(
            data[0], data[2],
            data[1], data[3]);
        }        

        static matrix<T,2> const Identity()
        {
            return matrix<T,2>(1.0,0.0,0.0,1.0);
        }

        static matrix<T,2> const Zero()
        {
            return matrix<T,2>(0.0,0.0,0.0,0.0);
        }

        static inline matrix<T,2> const Vecs2Mat(vec<T,2> const& a, vec<T,2> const& b)
        {
            return matrix<T,2>(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]);
        }

        inline matrix<T,2> normalized() const
        {
            T const det_ = determinant();
            return matrix<T,2>(data[0]/det_, data[1]/det_, 
                               data[2]/det_, data[3]/det_);
        }
        
        /* Standard operators */
        inline matrix<T,2> operator + (T const a) const
        {
            return matrix<T,2>(data[0]+a, data[1]+a,
                               data[2]+a, data[3]+a);
        }

        inline matrix<T,2> operator + (matrix<T,2> const& a) const
        {
            return matrix<T,2>(data[0]+a[0], data[1]+a[1],
                               data[2]+a[2], data[3]+a[3]);
        }

        inline matrix<T,2> operator - (T const a) const
        {
            return matrix<T,2>(data[0]-a, data[1]-a,
                               data[2]-a, data[3]-a);
        }

        inline matrix<T,2> operator - (matrix<T,2> const& a) const
        {
            return matrix<T,2>(data[0]-a[0], data[1]-a[1],
                               data[2]-a[2], data[3]-a[3]);
        }

        inline matrix<T,2> operator * (T const& a) const
        {
            return matrix<T,2>(data[0]*a, data[1]*a,
                               data[2]*a, data[3]*a);
        }

        inline vec<T,2> operator * (vec<T,2> const& a) const
        {
            return  vec<T,2>(data[0]*a[0] + data[1]*a[1],
                             data[2]*a[0] + data[3]*a[1]);
        }

        inline matrix<T,2> operator * (matrix<T,2> const& a) const
        {
            matrix<T,2> b;

            b(0,0) = data[0]*a[0] + data[1]*a[2];
            b(0,1) = data[0]*a[1] + data[1]*a[3];

            b(1,0) = data[2]*a[0] + data[3]*a[2];
            b(1,1) = data[2]*a[1] + data[3]*a[3];

            return b;
        }

        inline matrix<T,2> operator / (T const& a) const
        {
            return matrix<T,2>(data[0]/a, data[1]/a,
                               data[2]/a, data[3]/a);
        }

        inline matrix<T,2>& operator += (T const& a)
        {
            data[0] += a;
            data[1] += a;
            data[2] += a;
            data[3] += a;
            return *this;
        }

        inline matrix<T,2>& operator += (matrix<T,2> const& a)
        {
            data[0] += a[0];
            data[1] += a[1];
            data[2] += a[2];
            data[3] += a[3];
            return *this;
        }

        inline matrix<T,2>& operator -= (T const& a)
        {
            data[0] -= a;
            data[1] -= a;
            data[2] -= a;
            data[3] -= a;
            return *this;
        }

        inline matrix<T,2>& operator -= (matrix<T,2> const& a)
        {
            data[0] -= a[0];
            data[1] -= a[1];
            data[2] -= a[2];
            data[3] -= a[3];
            return *this;
        }

        inline matrix<T,2>& operator *= (T const& a)
        {
            data[0] *= a;
            data[1] *= a;
            data[2] *= a;
            data[3] *= a;
            return *this;
        }

        inline matrix<T,2>& operator *= (matrix<T,2> const& a) 
        {
            matrix<T,2> b;

            b(0,0) = data[0]*a[0] + data[1]*a[2];
            b(0,1) = data[0]*a[1] + data[1]*a[3];

            b(1,0) = data[2]*a[0] + data[3]*a[2];
            b(1,1) = data[2]*a[1] + data[3]*a[3];

            *this = b;
            return *this;
        }

        inline matrix<T,2>& operator /= (T const& a)
        {
            data[0] /= a;
            data[1] /= a;
            data[2] /= a;
            data[3] /= a;
            return *this;
        }

        T& operator()(int const row, int const col) { return data[row*2 + col]; }

        T const& operator () (int const row, int const col) const { return data[row*2 + col]; }

        size_t rows(){return 2;}
        size_t cols(){return 2;}
        size_t size(){return 4;}

    private:
        T operator [](size_t const index)  const
        {
            return data[index];
        }

        T data[4];
};

template <typename T>
class matrix<T,3>
{
    public:
        matrix(void) {}

        matrix<T,3>(T const a): data{a,a,a,a,a,a,a,a,a}{}


        matrix<T,3>(T const a, T const b, T const c,
               T const d, T const e, T const f,
               T const g, T const h, T const i):
             data{a,b,c,d,e,f,g,h,i} {}

        matrix<T,3>(T const a[9]): data{a} {}

        matrix<T,3>(matrix<T,3> const& a) { *this = a;}

        T const determinant() const
        {
            return 
            (data[0]*(data[4]*data[8] - data[7]*data[5]) - 
            data[1]*(data[3]*data[8] - data[6]*data[5]) +
            data[2]*(data[3]*data[7] - data[6]*data[4]));
        }

        matrix<T,3> transpose()
        {
            return matrix<T,3>(
            data[0], data[3], data[6],
            data[1], data[4], data[7],
            data[2], data[5], data[8]);
        }

        static matrix<T,3> const Identity()
        {
            return matrix<T,3>(1.0,0.0,0.0,
                               0.0,1.0,0.0,
                               0.0,0.0,1.0);
        }

        static matrix<T,3> const Zero()
        {
            return matrix<T,3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        }

        static matrix<T,3> const Vecs2Mat(vec<T,3> const& a, vec<T,3> const& b)
        {
            return matrix<T,3>(a[0]*b[0],a[0]*b[1],a[0]*b[2],
                            a[1]*b[0],a[1]*b[1],a[1]*b[2],
                            a[2]*b[0],a[2]*b[1],a[2]*b[2]);
        }

        inline matrix<T,3> normalized() const
        {
            T const det_ = determinant();
            return matrix<T,3>(data[0]/det_,data[1]/det_,data[2]/det_,
                               data[3]/det_,data[4]/det_,data[5]/det_,
                               data[6]/det_,data[7]/det_,data[8]/det_);
        }

        /* Standard operators */
        inline matrix<T,3> operator + (T const a) const
        {
            return matrix<T,3>(data[0]+a, data[1]+a, data[2]+a,
                               data[3]+a, data[4]+a, data[5]+a,
                               data[6]+a, data[7]+a, data[8]+a);
        }

        inline matrix<T,3> operator + (matrix<T,3> const& a) const
        {
            return matrix<T,3>(a[0]+data[0], a[1]+data[1], a[2]+data[2],
                               a[3]+data[3], a[4]+data[4], a[5]+data[5],
                               a[6]+data[6], a[7]+data[7], a[8]+data[8]);
        }

        inline matrix<T,3> operator - (T const a) const
        {
            return matrix<T,3>(data[0]-a, data[1]-a, data[2]-a,
                               data[3]-a, data[4]-a, data[5]-a,
                               data[6]-a, data[7]-a, data[8]-a);
        }

        inline matrix<T,3> operator - (matrix<T,3> const& a) const
        {
            return matrix<T,3>(data[0]-a[0], data[1]-a[1], data[2]-a[2],
                               data[3]-a[3], data[4]-a[4], data[5]-a[5],
                               data[6]-a[6], data[7]-a[7], data[8]-a[8]);
        }

        inline matrix<T,3> operator * (T const& a) const
        {
            return matrix<T,3>(data[0]*a, data[1]*a, data[2]*a,
                               data[3]*a, data[4]*a, data[5]*a,
                               data[6]*a, data[7]*a, data[8]*a);
        }

        inline vec<T,3> operator*(vec<T,3> const& a) const
        {
            return  vec<T,3>(data[0]*a[0] + data[1]*a[1] + data[2]*a[2],
                             data[3]*a[0] + data[4]*a[1] + data[5]*a[2],
                             data[6]*a[0] + data[7]*a[1] + data[8]*a[2]);
        }

        inline matrix<T,3> operator*(matrix<T,3> const& a) const
        {
            matrix<T,3> b;

            b(0,0) = data[0]*a[0] + data[1]*a[3] + data[2]*a[6];
            b(0,1) = data[0]*a[1] + data[1]*a[4] + data[2]*a[7];
            b(0,2) = data[0]*a[2] + data[1]*a[5] + data[2]*a[8];

            b(1,0) = data[3]*a[0] + data[4]*a[3] + data[5]*a[6];
            b(1,1) = data[3]*a[1] + data[4]*a[4] + data[5]*a[7];
            b(1,2) = data[3]*a[2] + data[4]*a[5] + data[5]*a[8];

            b(2,0) = data[6]*a[0] + data[7]*a[3] + data[8]*a[6];
            b(2,1) = data[6]*a[1] + data[7]*a[4] + data[8]*a[7];
            b(2,2) = data[6]*a[2] + data[7]*a[5] + data[8]*a[8];

            return b;
        }

        inline matrix<T,3> operator / (T const& a) const
        {
            return matrix<T,3>(data[0]/a, data[1]/a, data[2]/a,
                               data[3]/a, data[4]/a, data[5]/a,
                               data[6]/a, data[7]/a, data[8]/a);
        }

        inline matrix<T,3>& operator += (T const& a)
        {
            data[0] += a;
            data[1] += a;
            data[2] += a;
            data[3] += a;
            data[4] += a;
            data[5] += a;
            data[6] += a;
            data[7] += a;
            data[8] += a;
            return *this;
        }

        inline matrix<T,3>& operator += (matrix<T,3> const& a)
        {
            data[0] += a[0];
            data[1] += a[1];
            data[2] += a[2];
            data[3] += a[3];
            data[4] += a[4];
            data[5] += a[5];
            data[6] += a[6];
            data[7] += a[7];
            data[8] += a[8];
            return *this;
        }

        inline matrix<T,3>& operator -= (T const& a)
        {
            data[0] -= a;
            data[1] -= a;
            data[2] -= a;
            data[3] -= a;
            data[4] -= a;
            data[5] -= a;
            data[6] -= a;
            data[7] -= a;
            data[8] -= a;
            return *this;
        }

        inline matrix<T,3>& operator -= (matrix<T,3> const& a)
        {
            data[0] -= a[0];
            data[1] -= a[1];
            data[2] -= a[2];
            data[3] -= a[3];
            data[4] -= a[4];
            data[5] -= a[5];
            data[6] -= a[6];
            data[7] -= a[7];
            data[8] -= a[8];
            return *this;
        }

        inline matrix<T,3>& operator *= (T const& a)
        {
            data[0] *= a;
            data[1] *= a;
            data[2] *= a;
            data[3] *= a;
            data[4] *= a;
            data[5] *= a;
            data[6] *= a;
            data[7] *= a;
            data[8] *= a;
            return *this;
        }

        inline matrix<T,3>& operator /= (T const& a)
        {
            data[0] /= a;
            data[1] /= a;
            data[2] /= a;
            data[3] /= a;
            data[4] /= a;
            data[5] /= a;
            data[6] /= a;
            data[7] /= a;
            data[8] /= a;
            return *this;
        }

        size_t rows(){return 3;}
        size_t cols(){return 3;}
        size_t size(){return 9;}

        T& operator()(int const row, int const col) { return data[row*3 + col]; }

        T const& operator () (int const row, int const col) const { return data[row*3 + col]; }

    private:
        T operator [](size_t const index)  const
        {
            return data[index];
        }

        T data[9];
};


template <typename T>
class matrix<T,4>
{
    public:
        matrix<T,4>(){}

        matrix<T,4>(T const a): data{a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a} {}

        matrix<T,4>(T const a, T const b, T const c, T const d,
               T const e, T const f, T const g, T const h,
               T const i, T const j, T const k, T const l,
               T const m, T const n, T const o, T const p):
             data{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p} {}

        matrix<T,4>(T const a[16]): data{a} {}


        matrix<T,4>(matrix<T,4> const& a) { *this = a;}
        // {
        //     data[0] = a[0];   data[1] = a[1];   data[2] = a[2];   data[3] = a[3]; 
        //     data[4] = a[4];   data[5] = a[5];   data[6] = a[6];   data[7] = a[7];
        //     data[8] = a[8];   data[9] = a[9];   data[10] = a[10]; data[11] = a[11];
        //     data[12] = a[12]; data[13] = a[13]; data[14] = a[14]; data[15] = a[15];
        // }

        void row(int const index, vec<T,4> const& a)
        {
            data[index*4+0] = a[0];
            data[index*4+1] = a[1];
            data[index*4+2] = a[2];
            data[index*4+3] = a[3];
        }

        matrix<T,4> transpose() const
        {
            return matrix<T,4>(
            data[0], data[4], data[8], data[12],
            data[1], data[5], data[9], data[13],
            data[2], data[6], data[10], data[14],
            data[3], data[7], data[11], data[15]);
        }

        static matrix<T,4> const Identity()
        {
            return matrix<T,4>(1.0,0.0,0.0,0.0,
                               0.0,1.0,0.0,0.0,
                               0.0,0.0,1.0,0.0,
                               0.0,0.0,0.0,1.0);
        }

        static matrix<T,4> const Zero()
        {
            return matrix<T,4>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                               0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        }
        
        inline static matrix<T,4> const Vecs2Mat(vec<T,4> const& a, vec<T,4> const& b)
        {
            return matrix<T,4>(a[0]*b[0],a[0]*b[1],a[0]*b[2],a[0]*b[3],
                            a[1]*b[0],a[1]*b[1],a[1]*b[2],a[1]*b[3],
                            a[2]*b[0],a[2]*b[1],a[2]*b[2],a[2]*b[3],
                            a[3]*b[0],a[3]*b[1],a[3]*b[2],a[3]*b[3]);
        }

        inline matrix<T,4> normalized() const
        {
            T const det_ = determinant();
            return matrix<T,4>(data[0]/det_,data[1]/det_,data[2]/det_,data[3]/det_,
                               data[4]/det_,data[5]/det_,data[6]/det_,data[7]/det_,
                               data[8]/det_,data[9]/det_,data[10]/det_,data[11]/det_,
                               data[12]/det_,data[13]/det_,data[14]/det_,data[15]/det_);
        }

        // vec4& row(int const index)
        // {
        //     return vec4(*data[index*4],*data[index*4+1],*data[index*4+2],*data[index*4+3]);
        // }

        inline T const determinant() const
        {
            matrix<T,3> a,b,c,d;
            a = matrix<T,3>(data[5], data[6], data[7],
                     data[9], data[10],data[11],
                     data[13],data[14],data[15]);

            b = matrix<T,3>(data[4], data[6], data[7],
                     data[8], data[10],data[11],
                     data[12],data[14],data[15]);

            c = matrix<T,3>(data[4], data[5], data[7],
                     data[8], data[9], data[11],
                     data[12],data[13],data[15]);

            d = matrix<T,3>(data[4], data[5], data[6],
                     data[8], data[9], data[10],
                     data[12],data[13],data[14]);

            return 
            (data[0]*a.determinant() - data[1]*b.determinant() +
             data[2]*c.determinant() - data[3]*d.determinant());
        }

        // T determinant(void) const
        // {
        //     /*
        //     We need to create a temporary
        //     */
        //     matrix<T, 4> temp(*this);
        //     /*We convert the temporary to upper triangular form*/
        //     T det = T(1);
        //     for (unsigned c = 0; c < 4; ++c)
        //     {
        //         det = det*temp(c,c);
        //         for (unsigned r = c + 1; r < 4; ++r)
        //         {
        //             T ratio = temp(r, c) / temp(c, c);
        //             for (unsigned k = c; k < 4; k++)
        //             {
        //                 temp(r, k) = temp(r, k) - ratio * temp(c, k);
        //             }
        //         }
        //     }

        //     return det;
        // }

        inline matrix<T,4> operator*(matrix<T,4> const& a)
        {
            matrix<T,4> b;

            b(0,0) = data[0]*a[0] + data[1]*a[4] + data[2]*a[8]  + data[3]*a[12];
            b(0,1) = data[0]*a[1] + data[1]*a[5] + data[2]*a[9]  + data[3]*a[13];
            b(0,2) = data[0]*a[2] + data[1]*a[6] + data[2]*a[10] + data[3]*a[14];
            b(0,3) = data[0]*a[3] + data[1]*a[7] + data[2]*a[11] + data[3]*a[15];

            b(1,0) = data[4]*a[0] + data[5]*a[4] + data[6]*a[8]  + data[7]*a[12];
            b(1,1) = data[4]*a[1] + data[5]*a[5] + data[6]*a[9]  + data[7]*a[13];
            b(1,2) = data[4]*a[2] + data[5]*a[6] + data[6]*a[10] + data[7]*a[14];
            b(1,3) = data[4]*a[3] + data[5]*a[7] + data[6]*a[11] + data[7]*a[15];

            b(2,0) = data[8]*a[0] + data[9]*a[4] + data[10]*a[8]  + data[11]*a[12];
            b(2,1) = data[8]*a[1] + data[9]*a[5] + data[10]*a[9]  + data[11]*a[13];
            b(2,2) = data[8]*a[2] + data[9]*a[6] + data[10]*a[10] + data[11]*a[14];
            b(2,3) = data[8]*a[3] + data[9]*a[7] + data[10]*a[11] + data[11]*a[15];

            b(3,0) = data[12]*a[0] + data[13]*a[4] + data[14]*a[8]  + data[15]*a[12];
            b(3,1) = data[12]*a[1] + data[13]*a[5] + data[14]*a[9]  + data[15]*a[13];
            b(3,2) = data[12]*a[2] + data[13]*a[6] + data[14]*a[10] + data[15]*a[14];
            b(3,3) = data[12]*a[3] + data[13]*a[7] + data[14]*a[11] + data[15]*a[15];

            return b;
        }

        /* Standard operators */
        inline matrix<T,4> operator + (T const a) const
        {
            return matrix<T,4>(data[0]+a, data[1]+a, data[2]+a, data[3]+a,
                               data[4]+a, data[5]+a, data[6]+a, data[7]+a,
                               data[8]+a, data[9]+a, data[10]+a, data[11]+a,
                               data[12]+a, data[13]+a, data[14]+a, data[15]+a);
        }

        inline matrix<T,4> operator + (matrix<T,4> const& a) const
        {
            return matrix<T,4>(data[0]+a[0], data[1]+a[1], data[2]+a[2], data[3]+a[3],
                               data[4]+a[4], data[5]+a[5], data[6]+a[6], data[7]+a[7],
                               data[8]+a[8], data[9]+a[9], data[10]+a[10], data[11]+a[11],
                               data[12]+a[12], data[13]+a[13], data[14]+a[14], data[15]+a[15]);
        }

        inline matrix<T,4> operator - (T const a) const
        {
            return matrix<T,4>(data[0]-a, data[1]-a, data[2]-a, data[3]-a,
                               data[4]-a, data[5]-a, data[6]-a, data[7]-a,
                               data[8]-a, data[9]-a, data[10]-a, data[11]-a,
                               data[12]-a, data[13]-a, data[14]-a, data[15]-a);
        }

        inline matrix<T,4> operator - (matrix<T,4> const& a) const
        {
            return matrix<T,4>(data[0]-a[0], data[1]-a[1], data[2]-a[2], data[3]-a[3],
                               data[4]-a[4], data[5]-a[5], data[6]-a[6], data[7]-a[7],
                               data[8]-a[8], data[9]-a[9], data[10]-a[10], data[11]-a[11],
                               data[12]-a[12], data[13]-a[13], data[14]-a[14], data[15]-a[15]);
        }

        inline matrix<T,4> operator / (T const& a) const
        {
            return matrix<T,4>(data[0]/a, data[1]/a, data[2]/a, data[3]/a,
                               data[4]/a, data[5]/a, data[6]/a, data[7]/a,
                               data[8]/a, data[9]/a, data[10]/a, data[11]/a,
                               data[12]/a, data[13]/a, data[14]/a, data[15]/a);
        }

        inline matrix<T,4>& operator += (T const& a)
        {
            data[0] += a;
            data[1] += a;
            data[2] += a;
            data[3] += a;
            data[4] += a;
            data[5] += a;
            data[6] += a;
            data[7] += a;
            data[8] += a;
            data[9] += a;
            data[10] += a;
            data[11] += a;
            data[12] += a;
            data[13] += a;
            data[14] += a;
            data[15] += a;
            return *this;
        }

        inline matrix<T,4>& operator += (matrix<T,4> const& a)
        {
            data[0] += a[0];
            data[1] += a[1];
            data[2] += a[2];
            data[3] += a[3];
            data[4] += a[4];
            data[5] += a[5];
            data[6] += a[6];
            data[7] += a[7];
            data[8] += a[8];
            data[9] += a[9];
            data[10] += a[10];
            data[11] += a[11];
            data[12] += a[12];
            data[13] += a[13];
            data[14] += a[14];
            data[15] += a[15];
            return *this;
        }

        inline matrix<T,4>& operator -= (T const& a)
        {
            data[0] -= a;
            data[1] -= a;
            data[2] -= a;
            data[3] -= a;
            data[4] -= a;
            data[5] -= a;
            data[6] -= a;
            data[7] -= a;
            data[8] -= a;
            data[9] -= a;
            data[10] -= a;
            data[11] -= a;
            data[12] -= a;
            data[13] -= a;
            data[14] -= a;
            data[15] -= a;
            return *this;
        }

        inline matrix<T,4>& operator -= (matrix<T,4> const& a)
        {
            data[0] -= a[0];
            data[1] -= a[1];
            data[2] -= a[2];
            data[3] -= a[3];
            data[4] -= a[4];
            data[5] -= a[5];
            data[6] -= a[6];
            data[7] -= a[7];
            data[8] -= a[8];
            data[9] -= a[9];
            data[10] -= a[10];
            data[11] -= a[11];
            data[12] -= a[12];
            data[13] -= a[13];
            data[14] -= a[14];
            data[15] -= a[15];
            return *this;
        }

        inline matrix<T,4>& operator *= (T const& a)
        {
            data[0] *= a;
            data[1] *= a;
            data[2] *= a;
            data[3] *= a;
            data[4] *= a;
            data[5] *= a;
            data[6] *= a;
            data[7] *= a;
            data[8] *= a;
            data[9] *= a;
            data[10] *= a;
            data[11] *= a;
            data[12] *= a;
            data[13] *= a;
            data[14] *= a;
            data[15] *= a;
            return *this;
        }

        inline matrix<T,4>& operator /= (T const& a)
        {
            data[0] /= a;
            data[1] /= a;
            data[2] /= a;
            data[3] /= a;
            data[4] /= a;
            data[5] /= a;
            data[6] /= a;
            data[7] /= a;
            data[8] /= a;
            data[9] /= a;
            data[10] /= a;
            data[11] /= a;
            data[12] /= a;
            data[13] /= a;
            data[14] /= a;
            data[15] /= a;
            return *this;
        }


        unsigned const rows(){return 4;}
        unsigned const cols(){return 4;}
        unsigned const size(){return 16;}

        T& operator()(unsigned const row, unsigned const col) { return data[row*4 + col]; }

        T const& operator () (unsigned const row, unsigned const col) const { return data[row*4 + col]; }

    private:
        T operator [](unsigned const index)  const
        {
            return data[index];
        }

        T data[16];
};

/**************************************************************/
/*************** Full Pivoting LU Decomposition ***************/
/**************************************************************/
/* Traits: Rank revealing, providing a more robust check of   */
/* if a matrix is invertible. Eigenvalue revealing, which are */
/* also used. */
template<typename T, unsigned DIM> 
class FullPivLU
{
    static_assert(DIM > 0, "Dimension must be greater than 0");

    public:
    
      
    
        template<typename... Args>
        FullPivLU(Args&&... args) : data{T(args)...} {}

        T& operator()(int const row, int const col) { return data[row*DIM + col]; }

        T const& operator () (int const row, int const col) const { return data[row*DIM + col]; }

        size_t rows(){return DIM;}
        size_t cols(){return DIM;}
        size_t size(){return DIM*DIM;}

    private:
        T operator [](size_t const index)  const
        {
            return data[index];
        }

        T data[DIM*DIM];
        int p[DIM];
        int q[DIM];
        int nonzero_pivots, n_transp;
        T maxpivot, det_pq;
};

template<typename T, unsigned DIM>
FullPivLU<T,DIM> const operator*(T const& a, FullPivLU<T,DIM> const& b)
{
    return b*a;
}

template<typename T, unsigned DIM>
FullPivLU<T,DIM> const operator-(FullPivLU<T,DIM> const& b)
{
    return -b;
}

template <typename T>
class FullPivLU<T,2>
{
    public: 
        FullPivLU<T,2>(matrix<T,2> const& a)
        {
            /* Compute the pivot matrices */
            data[0] = a(0,0);
            data[1] = a(0,1);
            data[2] = a(1,0);
            data[3] = a(1,1);

            /*For  float, power is -24*/
            T MEPSILON = std::numeric_limits<T>::epsilon(); 

            n_transp = 0;
            nonzero_pivots = 2;
            maxpivot = 0.0;
            p[0] = 0; p[1] = 1;
            q[0] = 0; q[1] = 1;
            
            int ii,jj,kk;
            int pi,pj,tmp;
            T max;
            T ftmp;

            for (kk = 0; kk < 2; kk++)
            {
                pi=-1,pj=-1,max=0.0;
                //find pivot in submatrix *this->(k:n,k:n)
                for (ii = kk; ii < 2; ii++) 
                {
                    for (jj = kk; jj < 2; jj++) 
                    {
                        if (std::fabs(data[ii*2 + jj])>max)
                        {
                            max = std::fabs(data[ii*2 + jj]);
                            pi=ii;
                            pj=jj;
                        }
                    }
                }

                if(max < MEPSILON)
                {
                    nonzero_pivots = kk;
                    break;
                }

                if(max > maxpivot)
                {
                    maxpivot = max;
                }

                //Swap Row
                if(pi != kk)
                {
                    n_transp++;
                    tmp=p[kk];
                    p[kk]=p[pi];
                    p[pi]=tmp;
                    for (jj = 0; jj < 2; jj++)
                    {
                        ftmp = data[kk*2 + jj];
                        data[kk*2 + jj] = data[pi*2 + jj];
                        data[pi*2 + jj] = ftmp;
                    }
                }
                //Swap Col
                if(pj != kk)
                {
                    n_transp++;
                    tmp=q[kk];
                    q[kk]=q[pj];
                    q[pj]=tmp;
                    for (ii = 0; ii < 2; ii++)
                    {
                        ftmp = data[ii*2 + kk];
                        data[ii*2 + kk]= data[ii*2 + pj];
                        data[ii*2 + pj] = ftmp;
                    }
                }
                //END PIVOT
        
                //check pivot size and decompose
                if ((std::fabs(data[kk*2 + kk])> MEPSILON))
                {
                    for (ii = kk+1; ii < 2; ii++)
                    {
                        //Column normalisation
                        ftmp=data[ii*2 + kk] /= data[kk*2 + kk];
                        for (jj = kk+1; jj < 2; jj++)
                        {
                            //*this->(ik)**this->(kj) subtracted from lower right submatrix elements
                            data[ii*2 + jj] -= (ftmp * data[kk*2 + jj]);
                        }
                    }
                }
                //END DECOMPOSE
            }

            det_pq = (n_transp % 2) ? -1.0 : 1.0;
        }


        /* LU Decomposition functions */
        inline unsigned rank() const
        {
            T threshold = std::fabs(maxpivot) * std::numeric_limits<T>::epsilon() * 2.0;
            unsigned result = 0;
            for(int ii = 0; ii < nonzero_pivots; ++ii)
                result += (std::fabs(data[ii*2 + ii]) > threshold);
            return result;
        }

        inline bool isInvertible() const
        {
            return (rank() == 2);
        }

        inline matrix<T,2> const inverse() const
        {
            T det = determinant();
            return matrix<T,2>(data[3]/det, -data[1]/det, 
                               -data[2]/det, data[0]/det);
        }

        inline T determinant() const
        {
            return det_pq * data[0] * data[3];
        }

        inline T minEigenval() const
        {
            return std::min(data[0],data[3]);
        }

        inline T maxEigenval() const
        {
            return std::max(data[0],data[3]);
        }

        T& operator()(unsigned const row, unsigned const col) { return data[row*2 + col]; }

        T const& operator () (unsigned const row, unsigned const col) const { return data[row*2 + col]; }

        unsigned const rows() const {return 2;}
        unsigned const cols() const {return 2;}
        unsigned const size() const {return 4;}

    private:
        T operator [](unsigned const index)  const
        {
            return data[index];
        }

        T data[4];
        int p[2] = {0,0};
        int q[2] = {0,0};
        int nonzero_pivots, n_transp;
        T maxpivot, det_pq;
};


template <typename T>
class FullPivLU<T,3>
{
    public: 
        FullPivLU<T,3>(matrix<T,3> const& a)
        {
            /* Compute the pivot matrices */
            data[0] = a(0,0); data[1] = a(0,1); data[2] = a(0,2);
            data[3] = a(1,0); data[4] = a(1,1); data[5] = a(1,2);
            data[6] = a(2,0); data[7] = a(2,1); data[8] = a(2,2);

            /*For  float, power is -24*/
            T MEPSILON = std::numeric_limits<T>::epsilon(); 

            n_transp = 0;
            nonzero_pivots = 3;
            maxpivot = 0.0;
            
            p[0] = 0; p[1] = 1; p[2] = 2;
            q[0] = 0; q[1] = 1; q[2] = 2;

            int ii,jj,kk;
            int pi,pj,tmp;
            T max;
            T ftmp;

            for (kk = 0; kk < 3; kk++)
            {
                pi=-1,pj=-1,max=0.0;
                //find pivot in submatrix *this->(k:n,k:n)
                for (ii = kk; ii < 3; ii++) 
                {
                    for (jj = kk; jj < 3; jj++) 
                    {
                        if (std::fabs(data[ii*3 + jj])>max)
                        {
                            max = std::fabs(data[ii*3 + jj]);
                            pi=ii;
                            pj=jj;
                        }
                    }
                }

                if(max < MEPSILON)
                {
                    nonzero_pivots = kk;
                    break;
                }

                if(max > maxpivot)
                {
                    maxpivot = max;
                }

                //Swap Row
                if(pi != kk)
                {
                    n_transp++;
                    tmp=p[kk];
                    p[kk]=p[pi];
                    p[pi]=tmp;
                    for (jj = 0; jj < 3; jj++)
                    {
                        ftmp = data[kk*3 + jj];
                        data[kk*3 + jj] = data[pi*3 + jj];
                        data[pi*3 + jj] = ftmp;
                    }
                }
                //Swap Col
                if(pj != kk)
                {
                    n_transp++;
                    tmp=q[kk];
                    q[kk]=q[pj];
                    q[pj]=tmp;
                    for (ii = 0; ii < 3; ii++)
                    {
                        ftmp = data[ii*3 + kk];
                        data[ii*3 + kk]= data[ii*3 + pj];
                        data[ii*3 + pj] = ftmp;
                    }
                }
                //END PIVOT
        
                //check pivot size and decompose
                if ((std::fabs(data[kk*3 + kk])> MEPSILON))
                {
                    for (ii = kk+1; ii < 3; ii++)
                    {
                        //Column normalisation
                        ftmp=data[ii*3 + kk] /= data[kk*3 + kk];
                        for (jj = kk+1; jj < 3; jj++)
                        {
                            //*this->(ik)**this->(kj) subtracted from lower right submatrix elements
                            data[ii*3 + jj] -= (ftmp * data[kk*3 + jj]);
                        }
                    }
                }
                //END DECOMPOSE
            }

            det_pq = (n_transp % 2) ? -1.0 : 1.0;
        }

        /* LU Decomposition functions */
        inline unsigned rank() const
        {
            T threshold = abs(maxpivot) * std::numeric_limits<T>::epsilon() * 3.0;
            unsigned result = 0;
            for(unsigned ii = 0; ii < nonzero_pivots; ++ii)
                result += (std::abs(data[ii*3 + ii]) > threshold);
            return result;
        }

        inline bool isInvertible() const
        {
            return (rank() == 3);
        }

        // inline matrix<T,3> inverse() const
        // {
        //     T  det = determinant();
        //     return matrix<T,3>
        //     ((data[4]*data[8]-data[7]*data[5])/det,
        //     -(data[1]*data[8]-data[7]*data[2])/det,
        //      (data[1]*data[5]-data[2]*data[4])/det,

        //     -(data[3]*data[8]-data[6]*data[5])/det,
        //      (data[0]*data[8]-data[6]*data[2])/det,
        //     -(data[0]*data[5]-data[2]*data[4])/det,

        //      (data[3]*data[7]-data[4]*data[6])/det,
        //     -(data[0]*data[7]-data[1]*data[6])/det,
        //      (data[0]*data[4]-data[1]*data[3])/det);
        // }

        inline matrix<T,3> inverse() const
        {
            T const invdet = 1.0 / determinant();
            return matrix<T,3>(cofactor_3x3<0,0>(*this) * invdet,
                               cofactor_3x3<1,0>(*this) * invdet,
                               cofactor_3x3<2,0>(*this) * invdet,
                               cofactor_3x3<0,1>(*this) * invdet,
                               cofactor_3x3<1,1>(*this) * invdet,
                               cofactor_3x3<2,1>(*this) * invdet,
                               cofactor_3x3<0,2>(*this) * invdet,
                               cofactor_3x3<1,2>(*this) * invdet,
                               cofactor_3x3<2,2>(*this) * invdet);
        }

        inline T determinant() const
        {
            return det_pq * data[0] * data[4] * data[8];
        }

        /* The diagonal of the LU matrix is the eigenvalues */
        inline T minEigenval() const
        {
            return std::min(data[0],std::min(data[4],data[8]));
        }

        inline T maxEigenval() const
        {
            return std::max(data[0],std::max(data[4],data[8]));
        }

        T& operator()(unsigned const row, unsigned const col) { return data[row*3 + col]; }

        T const& operator () (unsigned const row, unsigned const col) const { return data[row*3 + col]; }

        unsigned const rows() const {return 3;}
        unsigned const cols() const {return 3;}
        unsigned const size() const {return 9;}

    private:

        template<int i, int j> 
        inline T const cofactor_3x3(FullPivLU<T,3> const& m) const
        {
            enum {
                i1 = (i+1) % 3,
                i2 = (i+2) % 3,
                j1 = (j+1) % 3,
                j2 = (j+2) % 3
            };

            return ((m(i1, j1) * m(i2, j2)) - (m(i1, j2) * m(i2, j1)));
        }

        T operator [](unsigned const index)  const
        {
            return data[index];
        }

        T data[9];
        unsigned p[3] = {0};
        unsigned q[3] = {0};
        unsigned nonzero_pivots, n_transp;
        T maxpivot, det_pq;
};

template <typename T>
class FullPivLU<T,4>
{
    public: 
        FullPivLU<T,4>(matrix<T,4> const& a)
        {
            /* Compute the pivot matrices */
            data[0] = a(0,0); data[1] = a(0,1); data[2] = a(0,2); data[3] = a(0,3);
            data[4] = a(1,0); data[5] = a(1,1); data[6] = a(1,2); data[7] = a(1,3);
            data[8] = a(2,0); data[9] = a(2,1); data[10] = a(2,2); data[11] = a(2,3);
            data[12] = a(3,0); data[13] = a(3,1); data[14] = a(3,2); data[15] = a(3,3);

            /*For  float, power is -24*/
            T MEPSILON = std::numeric_limits<T>::epsilon(); 

            n_transp = 0;
            nonzero_pivots = 4;
            maxpivot = 0.0;

            p[0] = 0; p[1] = 1; p[2] = 2; p[3] = 3;
            q[0] = 0; q[1] = 1; q[2] = 2; q[3] = 3;
            
            int ii,jj,kk;
            int pi,pj,tmp;
            T max;
            T ftmp;

            for (kk = 0; kk < 4; kk++)
            {
                pi=-1,pj=-1,max=0.0;
                //find pivot in submatrix *this->(k:n,k:n)
                for (ii = kk; ii < 4; ii++) 
                {
                    for (jj = kk; jj < 4; jj++) 
                    {
                        if (std::fabs(data[ii*4 + jj])>max)
                        {
                            max = std::fabs(data[ii*4 + jj]);
                            pi=ii;
                            pj=jj;
                        }
                    }
                }

                if(max < MEPSILON)
                {
                    nonzero_pivots = kk;
                    break;
                }

                if(max > maxpivot)
                {
                    maxpivot = max;
                }

                //Swap Row
                if(pi != kk)
                {
                    n_transp++;
                    tmp=p[kk];
                    p[kk]=p[pi];
                    p[pi]=tmp;
                    for (jj = 0; jj < 4; jj++)
                    {
                        ftmp = data[kk*4 + jj];
                        data[kk*4 + jj] = data[pi*4 + jj];
                        data[pi*4 + jj] = ftmp;
                    }
                }
                //Swap Col
                if(pj != kk)
                {
                    n_transp++;
                    tmp=q[kk];
                    q[kk]=q[pj];
                    q[pj]=tmp;
                    for (ii = 0; ii < 4; ii++)
                    {
                        ftmp = data[ii*4 + kk];
                        data[ii*4 + kk]= data[ii*4 + pj];
                        data[ii*4 + pj] = ftmp;
                    }
                }
                //END PIVOT
        
                //check pivot size and decompose
                if ((std::fabs(data[kk*4 + kk])> MEPSILON))
                {
                    for (ii = kk+1; ii < 4; ii++)
                    {
                        //Column normalisation
                        ftmp=data[ii*4 + kk] /= data[kk*4 + kk];
                        for (jj = kk+1; jj < 4; jj++)
                        {
                            //*this->(ik)**this->(kj) subtracted from lower right submatrix elements
                            data[ii*4 + jj] -= (ftmp * data[kk*4 + jj]);
                        }
                    }
                }
                //END DECOMPOSE
            }

            det_pq = (n_transp % 2) ? -1 : 1;
        }


        /* LU Decomposition functions */
        inline unsigned rank() const
        {
            T threshold = abs(maxpivot) * std::numeric_limits<T>::epsilon() * 4.0;
            unsigned result = 0;
            for(unsigned ii = 0; ii < nonzero_pivots; ++ii)
                result += (std::abs(data[ii*4 + ii]) > threshold);
            return result;
        }

        inline bool isInvertible() const
        {
            return (rank() == 4);
        }

        inline matrix<T,4> inverse() const
        {
            T det = determinant();
            matrix<T,4> inv = matrix<T,4>::Zero();

            inv(0,0) =  cofactor_4x4<0,0>(*this);
            inv(1,0) = -cofactor_4x4<0,1>(*this);
            inv(2,0) =  cofactor_4x4<0,2>(*this);
            inv(3,0) = -cofactor_4x4<0,3>(*this);
            inv(0,2) =  cofactor_4x4<2,0>(*this);
            inv(1,2) = -cofactor_4x4<2,1>(*this);
            inv(2,2) =  cofactor_4x4<2,2>(*this);
            inv(3,2) = -cofactor_4x4<2,3>(*this);
            inv(0,1) = -cofactor_4x4<1,0>(*this);
            inv(1,1) =  cofactor_4x4<1,1>(*this);
            inv(2,1) = -cofactor_4x4<1,2>(*this);
            inv(3,1) =  cofactor_4x4<1,3>(*this);
            inv(0,3) = -cofactor_4x4<3,0>(*this);
            inv(1,3) =  cofactor_4x4<3,1>(*this);
            inv(2,3) = -cofactor_4x4<3,2>(*this);
            inv(3,3) =  cofactor_4x4<3,3>(*this);

            return inv/=det;
        }

        inline T determinant() const
        {
            return det_pq * data[0] * data[5] * data[10] * data[15];
        }

        /* The diagonal of the LU matrix is the eigenvalues */
        inline T minEigenval() const
        {
            return std::min(data[0],std::min(data[5],std::min(data[10],data[15])));
        }

        inline T maxEigenval() const
        {
            return std::max(data[0],std::max(data[5],std::max(data[10],data[15])));
        }

        T& operator()(unsigned const row, unsigned const col) { return data[row*4 + col]; }

        T const& operator () (unsigned const row, unsigned const col) const { return data[row*4 + col]; }

        unsigned const rows() const {return 4;}
        unsigned const cols() const {return 4;}
        unsigned const size() const {return 16;}

    private:

        inline T const general_det3_helper
        (FullPivLU<T,4> const& mat, int i1, int i2, int i3, int j1, int j2, int j3)
        {
            return mat(i1,j1) * (mat(i2,j2) * mat(i3,j3) - mat(i2,j3) * mat(i3,j2));
        }

        template<int i, int j> 
        inline T const cofactor_4x4(FullPivLU<T,4> const& matrix)
        {
        enum {
            i1 = (i+1) % 4,
            i2 = (i+2) % 4,
            i3 = (i+3) % 4,
            j1 = (j+1) % 4,
            j2 = (j+2) % 4,
            j3 = (j+3) % 4
        };
        return general_det3_helper(matrix, i1, i2, i3, j1, j2, j3)
            + general_det3_helper(matrix, i2, i3, i1, j1, j2, j3)
            + general_det3_helper(matrix, i3, i1, i2, j1, j2, j3);
        }

        T operator [](unsigned const index)  const
        {
            return data[index];
        }

        T data[16];
        unsigned p[4] = {0};
        unsigned q[4] = {0};
        unsigned nonzero_pivots, n_transp;
        T maxpivot, det_pq;
};


#endif