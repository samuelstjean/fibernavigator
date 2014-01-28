#include "FMRIWindow.h"

#include "MainFrame.h"
#include "SceneManager.h"
#include "SelectionBox.h"
#include "SelectionEllipsoid.h"
#include "../main.h"
#include "../dataset/Anatomy.h"
#include "../dataset/Fibers.h"
#include "../dataset/RestingStateNetwork.h"
#include "../dataset/RTFMRIHelper.h"
#include "../misc/IsoSurface/CIsoSurface.h"
#include "../misc/IsoSurface/TriangleMesh.h"

#include <wx/checkbox.h>
#include <wx/grid.h>
#include <wx/tglbtn.h>
#include <wx/treectrl.h>


IMPLEMENT_DYNAMIC_CLASS( FMRIWindow, wxScrolledWindow )

BEGIN_EVENT_TABLE( FMRIWindow, wxScrolledWindow )
END_EVENT_TABLE()


FMRIWindow::FMRIWindow( wxWindow *pParent, MainFrame *pMf, wxWindowID id, const wxPoint &pos, const wxSize &size)
:   wxScrolledWindow( pParent, id, pos, size, wxBORDER_NONE, _T("fMRI resting-state networks") ),
    m_pMainFrame( pMf ),
	showRawData( true )
{
    SetBackgroundColour( *wxLIGHT_GREY );
    m_pFMRISizer = new wxBoxSizer( wxVERTICAL );
    SetSizer( m_pFMRISizer );
    SetAutoLayout( true );

    m_pBtnSelectFMRI = new wxButton( this, wxID_ANY,wxT("Load time-courses"), wxDefaultPosition, wxSize(230, -1) );
	pMf->Connect( m_pBtnSelectFMRI->GetId(), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MainFrame::onLoadAsRestingState) );
    m_pBtnSelectFMRI->SetBackgroundColour(wxColour( 255, 147, 147 ));

    m_pBtnSelectClusters = new wxButton( this, wxID_ANY,wxT("Clusters not selected"), wxDefaultPosition, wxSize(230, -1) );
    Connect( m_pBtnSelectClusters->GetId(), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(FMRIWindow::OnSelectClusters) );
    m_pBtnSelectClusters->SetBackgroundColour(wxColour( 255, 147, 147 ));
	m_pBtnSelectClusters->Enable(false);

	m_pBtnStart = new wxToggleButton( this, wxID_ANY,wxT("Start correlation"), wxDefaultPosition, wxSize(230, 50) );
    Connect( m_pBtnStart->GetId(), wxEVT_COMMAND_TOGGLEBUTTON_CLICKED, wxCommandEventHandler(FMRIWindow::OnStartRTFMRI) );
    m_pBtnStart->Enable(false);

	wxBoxSizer *pBoxRow1 = new wxBoxSizer( wxVERTICAL );
	pBoxRow1->Add( m_pBtnSelectFMRI, 0, wxALIGN_CENTER | wxALL, 1 );
    pBoxRow1->Add( m_pBtnSelectClusters, 0, wxALIGN_CENTER, 1 );
	pBoxRow1->Add( m_pBtnStart, 0, wxALIGN_CENTER | wxALL, 1 );
	m_pFMRISizer->Add( pBoxRow1, 0, wxFIXED_MINSIZE | wxALL, 2 );

	m_pTextCorrThreshold = new wxStaticText( this, wxID_ANY, wxT("Z-Threshold"), wxDefaultPosition, wxSize(70, -1), wxALIGN_CENTER );
	m_pSliderCorrThreshold = new MySlider( this, wxID_ANY, 0, 0, 500, wxDefaultPosition, wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
	m_pSliderCorrThreshold->SetValue( 300 );
	Connect( m_pSliderCorrThreshold->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(FMRIWindow::OnSliderCorrThreshMoved) );
    m_pTxtCorrThreshBox = new wxTextCtrl( this, wxID_ANY, wxT("3.0"), wxDefaultPosition, wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	wxBoxSizer *pBoxRow5 = new wxBoxSizer( wxHORIZONTAL );
    pBoxRow5->Add( m_pTextCorrThreshold, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxRow5->Add( m_pSliderCorrThreshold,   0, wxALIGN_LEFT | wxEXPAND | wxALL, 1);
	pBoxRow5->Add( m_pTxtCorrThreshBox,   0, wxALIGN_LEFT | wxALL, 1);
	m_pFMRISizer->Add( pBoxRow5, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );

	m_pTextSizeP = new wxStaticText( this, wxID_ANY, wxT("Point size"), wxDefaultPosition, wxSize(70, -1), wxALIGN_CENTER );
	m_pSliderSizeP = new MySlider( this, wxID_ANY, 0, 1, 100, wxDefaultPosition, wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
	m_pSliderSizeP->SetValue( 10 );
	Connect( m_pSliderSizeP->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(FMRIWindow::OnSliderSizePMoved) );
    m_pTxtSizePBox = new wxTextCtrl( this, wxID_ANY, wxT("10.0"), wxDefaultPosition, wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	wxBoxSizer *pBoxRow7 = new wxBoxSizer( wxHORIZONTAL );
    pBoxRow7->Add( m_pTextSizeP, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxRow7->Add( m_pSliderSizeP,   0, wxALIGN_LEFT | wxEXPAND | wxALL, 1);
	pBoxRow7->Add( m_pTxtSizePBox,   0, wxALIGN_LEFT | wxALL, 1);
	m_pFMRISizer->Add( pBoxRow7, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );

	m_pTextAlpha = new wxStaticText( this, wxID_ANY, wxT("Alpha blend"), wxDefaultPosition, wxSize(70, -1), wxALIGN_CENTER );
	m_pSliderAlpha = new MySlider( this, wxID_ANY, 0, 0, 100, wxDefaultPosition, wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
	m_pSliderAlpha->SetValue( 50 );
	Connect( m_pSliderAlpha->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(FMRIWindow::OnSliderAlphaMoved) );
    m_pTxtAlphaBox = new wxTextCtrl( this, wxID_ANY, wxT("0.8"), wxDefaultPosition, wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	wxBoxSizer *pBoxRow8 = new wxBoxSizer( wxHORIZONTAL );
    pBoxRow8->Add( m_pTextAlpha, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxRow8->Add( m_pSliderAlpha,   0, wxALIGN_LEFT | wxEXPAND | wxALL, 1);
	pBoxRow8->Add( m_pTxtAlphaBox,   0, wxALIGN_LEFT | wxALL, 1);
	m_pFMRISizer->Add( pBoxRow8, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );

	m_pBtnConvertFMRI = new wxButton( this, wxID_ANY,wxT("Convert to Overlay"), wxDefaultPosition, wxSize(230, -1) );
	Connect( m_pBtnConvertFMRI->GetId(), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(FMRIWindow::onConvertRestingState) );

	wxBoxSizer *pBoxRow9 = new wxBoxSizer( wxHORIZONTAL );
	pBoxRow9->Add( m_pBtnConvertFMRI,   0, wxALIGN_LEFT | wxALL, 1);
	m_pFMRISizer->Add( pBoxRow9, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );

}

void FMRIWindow::OnSize( wxSizeEvent &WXUNUSED(event) )
{
	
}

void FMRIWindow::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
    wxPaintDC dc( this );
}

wxSizer* FMRIWindow::getWindowSizer()
{
    return m_pFMRISizer;
}

void FMRIWindow::OnSelectClusters( wxCommandEvent& WXUNUSED(event) )
{
	//Select map for threshold seeding
    long item = m_pMainFrame->getCurrentListIndex();
	Anatomy* pMap = (Anatomy*)DatasetManager::getInstance()->getDataset (MyApp::frame->m_pListCtrl->GetItem( item )); 

	if( pMap != NULL && pMap->getBands() == 1 )
    {
		m_pBtnSelectClusters->SetLabel( pMap->getName() );
        m_pBtnSelectClusters->SetBackgroundColour(wxNullColour);
        m_pMainFrame->m_pRestingStateNetwork->setClusters( (Anatomy *)DatasetManager::getInstance()->getDataset( m_pMainFrame->m_pListCtrl->GetItem( item ) ) );
	}
	
	m_pBtnStart->Enable( true ); 
	m_pBtnStart->SetBackgroundColour(wxColour( 147, 255, 239 ));
}

void FMRIWindow::SetSelectButton()
{
	m_pBtnSelectFMRI->SetLabel( wxT("Time-courses OK"));
    m_pBtnSelectFMRI->SetBackgroundColour(wxNullColour);

	m_pBtnSelectClusters->Enable(true);
}

void FMRIWindow::OnSliderCorrThreshMoved(  wxCommandEvent& WXUNUSED(event) )
{
	float sliderValue = m_pSliderCorrThreshold->GetValue() / 100.0f;
    m_pTxtCorrThreshBox->SetValue( wxString::Format( wxT( "%.2f"), sliderValue ) );
	m_pMainFrame->m_pRestingStateNetwork->SetCorrThreshold( sliderValue );
	RTFMRIHelper::getInstance()->setRTFMRIDirty( true );
}

void FMRIWindow::OnSliderSizePMoved(  wxCommandEvent& WXUNUSED(event) )
{
	float sliderValue = m_pSliderSizeP->GetValue();
	m_pTxtSizePBox->SetValue( wxString::Format( wxT( "%.1f"), sliderValue ) );
	m_pMainFrame->m_pRestingStateNetwork->SetSizePSliderValue( sliderValue );
	RTFMRIHelper::getInstance()->setRTFMRIDirty( true );
}


void FMRIWindow::OnSliderAlphaMoved(  wxCommandEvent& WXUNUSED(event) )
{
	float sliderValue = m_pSliderAlpha->GetValue() / 100.0f;
	m_pTxtAlphaBox->SetValue( wxString::Format( wxT( "%.2f"), sliderValue ) );
	m_pMainFrame->m_pRestingStateNetwork->SetAlphaSliderValue( sliderValue );
	RTFMRIHelper::getInstance()->setRTFMRIDirty( true );
}

void FMRIWindow::onConvertRestingState( wxCommandEvent& WXUNUSED(event) )
{
	//Convert to anat
	m_pMainFrame->m_pRestingStateNetwork->SetExport( true );
	m_pMainFrame->m_pRestingStateNetwork->seedBased();
	m_pMainFrame->m_pRestingStateNetwork->SetExport( false );

	
	int indx = DatasetManager::getInstance()->createAnatomy( m_pMainFrame->m_pRestingStateNetwork->getZscores(), OVERLAY );
    
	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( indx );
    pNewAnatomy->setShowFS(false);

    pNewAnatomy->setType(OVERLAY);
    pNewAnatomy->setDataType(16);
    pNewAnatomy->setName( wxT("Z-score map") );
    MyApp::frame->m_pListCtrl->InsertItem( indx );

	RTFMRIHelper::getInstance()->setRTFMRIReady(false);

	m_pMainFrame->m_pRestingStateNetwork->clear3DPoints();
	RTFMRIHelper::getInstance()->setRTFMRIDirty( false );
    m_pBtnStart->SetLabel(wxT("Start correlation"));
    m_pBtnStart->SetValue(false);

}

void FMRIWindow::OnStartRTFMRI( wxCommandEvent& WXUNUSED(event) )
{
	RTFMRIHelper::getInstance()->toggleRTFMRIReady();
    RTFMRIHelper::getInstance()->setRTFMRIDirty( true );

    if( !RTFMRIHelper::getInstance()->isRTFMRIReady() )
    {
		m_pMainFrame->m_pRestingStateNetwork->clear3DPoints();
        RTFMRIHelper::getInstance()->setRTFMRIDirty( false );
        m_pBtnStart->SetLabel(wxT("Start correlation"));
    }
    else
    {
        m_pBtnStart->SetLabel(wxT("Stop correlation"));
	}
}

