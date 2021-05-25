#include <iostream>
#include <vector>
#include <WeakPtrFactories.h>
#include <StaticID.h>

using namespace model;

struct MyContainer;

struct MyClass
{
    
    typedef size_t KeyType;
    const KeyType key;
    
    MyClass(MyContainer* const cont,const KeyType& key_in) :
    key(key_in)
    {
        
    }
    
};

struct MyContainer : public WeakPtrFactory<MyContainer,MyClass>
{
    
};





int main()
{

    MyContainer c;
    
//    std::vector<std::shared_ptr<MyClass>> v;
    
    auto m4(c.create(2));
//
//
    {
//        auto m0(c.create(0));
    }
//
//    auto m3(c.create(3));

    
//    for(int k=1;k<300;++k)
//    {
//        v.push_back(c.create(k));
//    }
//    auto m1(c.create(1));
//    auto m0(c.create(0));
    {
        auto m1(c.create(1));
                auto m2(c.create(3));
    }
  //  auto m2(c.create(2));

    
//    for(typename MyContainer::iterator pair=c.begin();pair!=c.end();++pair)
//    {
//        std::cout<<pair->first<<" "<<pair->second.lock()->key<<" "<<c.keyID(pair->first)<<std::endl;
////        std::cin.get();
//    }
    
    std::cout<<"c.size="<<c.size()<<std::endl;
    
    for(const auto& pair : c)
    {
        std::cout<<pair.first<<" "<<pair.second.lock()->key<<" "<<c.keyID(pair.first)<<std::endl;
        //        std::cin.get();
    }
    
//    c.keyID(5);
    
    return 0;
}
