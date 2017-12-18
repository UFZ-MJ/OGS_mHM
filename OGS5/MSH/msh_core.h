/**************************************************************************
MSHLib - Object:
Task:
Programing:
08/2005 WW/OK Encapsulation from rf_ele_msh
last modified
**************************************************************************/
#ifndef msh_core_INC
#define msh_core_INC

#include <string>
#include <iostream>

//------------------------------------------------------------------------
namespace FiniteElement
{
   class CElement;
   class CFiniteElementStd;
   class CFiniteElementVec;
   class ElementMatrix;
   class ElementMatrix_DM;
}


namespace Mesh_Group
{

   //------------------------------------------------------------------------
   class CCore
   {
      public:
         // Methods
         CCore(size_t id);
         virtual ~CCore() {}

         // Get members
         size_t GetIndex() const { return index; }

         bool GetMark() const
         {
            return mark;
         }

         bool GetOrder() const { return quadratic; };

         bool Interior() const
         {
            if (boundary_type == 'I')
               return true;
            else
               return false;
         }

         // Set members
         void SetBoundaryType(char type)
         {
            boundary_type = type;
         }

         char GetBoundaryType() const
         {
            return boundary_type;
         }                                        // 18.02.2009. WW

         void SetOrder(bool order)
         {
            quadratic = order;
         }

         void SetMark(bool state)
         {
            mark = state;
         }

         void SetIndex(size_t lvalue)
         {
            index = lvalue;
         }                                        //OK

         // Output
         virtual void Write(std::ostream& os = std::cout) const
         {
            os << std::endl;
         }

      protected:
         // Properties
         size_t index;
         char boundary_type;
         bool mark;                               // e.g. marked to be refined or active
         bool quadratic;                          // higher order

   };

}                                                 // namespace Mesh_Group
#endif
