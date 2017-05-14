/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_vtkActorWrapper_H_
#define model_vtkActorWrapper_H_

// http://stackoverflow.com/questions/20945017/how-to-subclass-vtkactor

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkObjectFactory.h>
#include <vtkRenderingCoreModule.h>
#include <vtkProperty.h>


namespace model
{
    
    
    
//    class T;
     template <typename T>
    class VTKRENDERINGCORE_EXPORT vtkActorWrapper : public vtkActor {
    public:
        vtkTypeMacro(vtkActorWrapper, vtkActor);
        
        static vtkActorWrapper *New();
        
        virtual void ReleaseGraphicsResources(vtkWindow *window) {
            this->Device->ReleaseGraphicsResources(window);
            this->Superclass::ReleaseGraphicsResources(window);
        }
        
        virtual int RenderOpaqueGeometry(vtkViewport *viewport){
            if ( ! this->Mapper ) {
                return 0;
            }
            if (!this->Property) {
                this->GetProperty();
            }
            if (this->GetIsOpaque()) {
                vtkRenderer *ren = static_cast<vtkRenderer *>(viewport);
                this->Render(ren);
                return 1;
            }
            return 0;
        }
        
        virtual int RenderTranslucentPolygonalGeometry(vtkViewport *viewport){
            if ( ! this->Mapper ) {
                return 0;
            }
            if (!this->Property) {
                this->GetProperty();
            }
            if (!this->GetIsOpaque()) {
                vtkRenderer *ren = static_cast<vtkRenderer *>(viewport);
                this->Render(ren);
                return 1;
            }
            return 0;
        }
        
        virtual void Render(vtkRenderer *ren){
            this->Property->Render(this, ren);
            this->Device->SetProperty (this->Property);
            this->Property->Render(this, ren);
            if (this->BackfaceProperty) {
                this->BackfaceProperty->BackfaceRender(this, ren);
                this->Device->SetBackfaceProperty(this->BackfaceProperty);
            }
            if (this->Texture) {
                this->Texture->Render(ren);
            }
            this->ComputeMatrix();
            this->Device->SetUserMatrix(this->Matrix);
            this->Device->Render(ren,this->Mapper);
        }
        
        void ShallowCopy(vtkProp *prop) {
            vtkActorWrapper *f = vtkActorWrapper::SafeDownCast(prop);
            this->vtkActor::ShallowCopy(prop);
        }
        
        //****************************************//
        //              my member
        //****************************************//
        T*   userData{nullptr};
        
    protected:
        vtkActor* Device;
        
        vtkActorWrapper() {
            this -> Device = vtkActor::New();
        }
        
        ~vtkActorWrapper() {
            this -> Device -> Delete();
        }
    private:
    };
    
    vtkStandardNewMacro(vtkActorWrapper)
    
    
} // namespace model
#endif







