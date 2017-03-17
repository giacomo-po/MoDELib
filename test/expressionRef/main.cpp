#include <iostream>

#include <utility> // for
#include <memory> // shared_ptr
//#include <model/FEM/TrialOperators/ExpressionRef.h>

//template <typename T>
//struct ExpressionRef
//{
//    
//    const T& ref;
//
//    ExpressionRef(const T& ) :
//    
//};


template <typename Derived>
struct ScalarExp //: std::enable_shared_from_this<ScalarExp<Derived>>
{

    ScalarExp(){}
    
    ScalarExp(const ScalarExp<Derived>& other)=delete;
    ScalarExp(ScalarExp<Derived>&& other)=default;

    
    /**********************************************************************/
    const Derived& derived() const
    {/*! A const reference to the Derived object
      */
        return *static_cast<const Derived*>(this);
    }
    
    /**********************************************************************/
    Derived& derived()
    {/*! A const reference to the Derived object
      */
        return *static_cast<Derived*>(this);
    }
    
    ~ScalarExp()
    {
        std::cout<<"destroying ScalarExp "<<typeid(derived()).name()<<std::endl;

    }
    
//    std::shared_ptr<Derived> shared()
//    {
//        return this->shared_from_this();
//    }
    
//    std::shared_ptr<const Derived> shared() const
//    {
//        return std::static_pointer_cast<const Derived>(this->shared_from_this());
//    }
    
//    std::shared_ptr<const Derived> shared() const
//    {
//        return std::static_pointer_cast<const Derived>(this->shared_from_this());
//    }
    
//    double val() const
//    {
//        return derived().val();
//    }

};

struct Dref
{
const double& v;

Dref(const double& v_in) :v(v_in){}

};

struct Scalar : public ScalarExp<Scalar>
//std::enable_shared_from_this<Scalar>
{

//    const double& v;
//    
//    Scalar(const double& v_in) :v(v_in){}

    
    const Dref& dref;
        Scalar(const Dref& v_in) :dref(v_in){}

    
    template<typename T>
    Scalar(const ScalarExp<T>& ex) :dref(ex.val()){}

    Scalar(const Scalar& other)=delete;

//    Scalar(const Scalar& other) : v(other.v)
//    {
//        std::cout<<"Scalar Copy"<<std::endl;
//    }
    
    Scalar(Scalar&& other)=default;

    
    
//    double& val()
//    {
//        return v;
//    }

    const double& val() const
    {
        return dref.v;
    }

//    std::shared_ptr<const Scalar>&& shared() const
//    {
//        return this->shared_from_this();
//    }
//    
//    std::shared_ptr<Scalar>&& shared()
//    {
//        return this->shared_from_this();
//    }
    
};

template <typename T>
struct ScalarProd : ScalarExp<ScalarProd<T>>
{

    const double& op1;
    const T& op2;

    ScalarProd(const double& d,const T& exp) :
    op1(d),
    op2(exp)
    {
    
    }
    
    double val() const
    {
        return op1*op2.val();
    }
    
};

template<typename T>
ScalarProd<T> operator*(const double& d,
                           const ScalarExp<T>& ex)
{
    std::cout<<"operator* 1"<<std::endl;
    return ScalarProd<T>(d,ex.derived());
}

template<typename T1,typename T2>
class ScalarSum : public ScalarExp<ScalarSum<T1,T2>>
{


    
public:

    
    
    
//    const model::internal::ExpressionRef<T1> op1;
//    const model::internal::ExpressionRef<T2> op2;
    
//
//    ScalarSum(const std::shared_ptr<T1>& sp1,
//              const std::shared_ptr<T2>& sp2):
//    /* copy construct to share ownership*/ op1(sp1),
//    /* copy construct to share ownership*/ op2(sp2)
//    {
//        std::cout<<"ScalarSum Constructor 1"<<std::endl;
//    }
//    
//    double val() const
//    {
//        return op1->val()+op2->val();
//    }
    
    const std::shared_ptr<T1> p1;
    const std::shared_ptr<T2> p2;
 
    const T1& op1;
    const T2& op2;

    //
    ScalarSum(const T1& t1,
              const T2& t2) :
    p1(nullptr),
    p2(nullptr),
    op1(t1),
    op2(t2)
    {
        std::cout<<"ScalarSum Constructor 1"<<std::endl;
        std::cout<<typeid(t1).name()<<std::endl;
        std::cout<<typeid(t2).name()<<std::endl;
    }
    
//    ScalarSum(T1&& t1,
//              const T2& t2) :
//    p1(nullptr),
//    p2(nullptr),
//    op1(t1),
//    op2(t2)
//    {
//        std::cout<<"ScalarSum Constructor 2"<<std::endl;
//        std::cout<<typeid(t1).name()<<std::endl;
//        std::cout<<typeid(t2).name()<<std::endl;
//    }
//    
//    ScalarSum(const T1& t1,
//            T2&& t2) :
//    p1(nullptr),
//    p2(nullptr),
//    op1(t1),
//    op2(t2)
//    {
//        std::cout<<"ScalarSum Constructor 3"<<std::endl;
//        std::cout<<typeid(t1).name()<<std::endl;
//        std::cout<<typeid(t2).name()<<std::endl;
//    }
//    
//    ScalarSum(T1&& t1,
//              T2&& t2) :
//    p1(nullptr),
//    p2(nullptr),
//    op1(t1),
//    op2(t2)
//    {
//        std::cout<<"ScalarSum Constructor 4"<<std::endl;
//        std::cout<<typeid(t1).name()<<std::endl;
//        std::cout<<typeid(t2).name()<<std::endl;
//    }
//
    double val() const
    {
        return op1.val()+op2.val();
    }
    

    
    
};

//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const ScalarExp<T1>& ex1,
//                           const ScalarExp<T2>& ex2)
//{
//    std::cout<<"operator 1"<<std::endl;
//    return ScalarSum<T1,T2>(ex1.shared(),ex2.shared());
//}

//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const std::shared_ptr<T1>& p1,
//                           const std::shared_ptr<T2>& p2)
//{
//    std::cout<<"operator 1"<<std::endl;
//    return ScalarSum<T1,T2>(p1,p2);
//}
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const T1& p1,
//                           const T2& p2)
//{
//    std::cout<<"operator 2"<<std::endl;
//    return ScalarSum<T1,T2>(p1,p2);
//}

//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const ScalarExp<T1>& p1,
//                           const ScalarExp<T2>& p2)
//{
//    std::cout<<"operator 2"<<std::endl;
//    return ScalarSum<T1,T2>(p1,p2);
//}


//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const ScalarExp<T1>& p1,
//                           const ScalarExp<T2>& p2)
//{
//    std::cout<<"operator 1"<<std::endl;
//    return ScalarSum<T1,T2>(std::make_shared<T1>(p1.derived()),
//                            std::make_shared<T2>(p2.derived()));
//}
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(ScalarExp<T1>&& p1,
//                           const ScalarExp<T2>& p2)
//{
//    std::cout<<"operator 2"<<std::endl;
//    return ScalarSum<T1,T2>(std::make_shared<T1>(std::move(p1.derived())),
//                            std::make_shared<T2>(p2.derived()));
//}
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const ScalarExp<T1>& p1,
//                        ScalarExp<T2>&& p2)
//{
//    std::cout<<"operator 3"<<std::endl;
//    return ScalarSum<T1,T2>(std::make_shared<T1>(p1.derived()),
//                            std::make_shared<T2>(std::move(p2.derived())));
//}
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(ScalarExp<T1>&& p1,
//                           ScalarExp<T2>&& p2)
//{
//    std::cout<<"operator 4"<<std::endl;
//    return ScalarSum<T1,T2>(std::make_shared<T1>(std::move(p1.derived())),
//                            std::make_shared<T2>(std::move(p2.derived())));
//}

template<typename T1,typename T2>
ScalarSum<T1,T2> operator+(const ScalarExp<T1>& ex1,
                           const ScalarExp<T2>& ex2)
{
    std::cout<<"operator 1"<<std::endl;
    return ScalarSum<T1,T2>(ex1.derived(),ex2.derived());
}

//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(ScalarExp<T1>&& ex1,
//                           const ScalarExp<T2>& ex2)
//{
//    std::cout<<"operator 2"<<std::endl;
//    return ScalarSum<T1,T2>(std::move(ex1.derived()),ex2.derived());
//}
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const ScalarExp<T1>& ex1,
//                        ScalarExp<T2>&& ex2)
//{
//    std::cout<<"operator 3"<<std::endl;
//    return ScalarSum<T1,T2>(ex1.derived(),std::move(ex2.derived()));
//}
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(ScalarExp<T1>&& ex1,
//                           ScalarExp<T2>&& ex2)
//{
//    std::cout<<"operator 4"<<std::endl;
//    return ScalarSum<T1,T2>(std::move(ex1.derived()),std::move(ex2.derived()));
//}

//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(ScalarExp<T1>&& ex1,
//                           const ScalarExp<T2>& ex2)
//{
//        std::cout<<"operator 2"<<std::endl;
//    return ScalarSum<T1,T2>(std::move(ex1.derived()),ex2.derived());
//}
//
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(const ScalarExp<T1>& ex1,
//                           ScalarExp<T2>&& ex2)
//{
//        std::cout<<"operator 3"<<std::endl;
//    return ScalarSum<T1,T2>(ex1.derived(),std::move(ex2.derived()));
//}
//
//template<typename T1,typename T2>
//ScalarSum<T1,T2> operator+(ScalarExp<T1>&& ex1,
//                           ScalarExp<T2>&& ex2)
//{
//        std::cout<<"operator 4"<<std::endl;
//    return ScalarSum<T1,T2>(std::move(ex1.derived()),std::move(ex2.derived()));
//}

//struct OverScalar
//{
//    const Scalar& s;
//
//    OverScalar(const Scalar& sin) : s(sin){}
//
//};

int main()
{
    
//    Scalar s1(1.0);
//    Scalar s2(2.0);

    double d=4.5;
//    auto s3=(Scalar(2.0)+Scalar(d));
//
//    std::cout<<"s3.op1.val()="<<s3.op1.val()<<std::endl;
//    std::cout<<"s3.op2.val()="<<s3.op2.val()<<std::endl;
//    std::cout<<"s3.val()="<<s3.val()<<std::endl;
    
//    OverScalar os(Scalar(9.0));
//    
//    std::cout<<"os.s.val()="<<os.s.val()<<std::endl;
//
//    auto scalarsum=
//    std::shared_ptr<Scalar>(new Scalar(1+3.4))
//    +std::shared_ptr<Scalar>(new Scalar(2))
//    ;
    
//    auto scalarsum=
//    Scalar(1+3.4)+
//    Scalar(2)
//    ;

//    ScalarSum<Scalar,Scalar> scalarsum( Scalar(1+3.4),Scalar(2));

//    const Scalar s2(3.0);
    
    
//    const std::shared_ptr<Scalar> sp(std::make_shared<Scalar>(s2));
    
 //   std::cout<<sp->val()<<std::endl;

    
    auto scalarsum = 2.0*Scalar(Dref(2+3))+Scalar(Dref(3.0))+3.5*Scalar(Dref(1.0));
//    std::cout<<"scalarsum.val()="<<scalarsum.op1.op1.val()<<std::endl;
//    std::cout<<"scalarsum.val()="<<scalarsum.op1.op2.val()<<std::endl;
    std::cout<<"scalarsum.op1.val()="<<scalarsum.op1.val()<<std::endl;
    std::cout<<"scalarsum.op2.val()="<<scalarsum.op2.val()<<std::endl;
    std::cout<<"scalarsum.val()="<<scalarsum.val()<<std::endl;
//    std::cout<<"&scalarsum.op1.val()="<<&scalarsum.op1<<std::endl;
//    std::cout<<"&scalarsum.op2.val()="<<&scalarsum.op2<<std::endl;

    //
//    int i=6;
//    //A b(i);
//    //A b(i);
//    
//    B b((A(i)));
//
//    
//    std::cout<<"a.i="<<b.a.i<<std::endl;
//    
////    const int& j(3+8);
////
////        std::cout<<"j="<<j<<std::endl;
    
    return 0;
}
